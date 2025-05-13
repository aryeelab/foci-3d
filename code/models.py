import pytorch_lightning as pl
import torch
import math
from torch import nn

class UNetLightning(pl.LightningModule):
    def __init__(
        self,
        in_channels: int = 1,
        base_filters: int = 32,
        lr: float = 1e-3,
        weight_decay: float = 1e-5
    ):
        super().__init__()
        self.save_hyperparameters()

        # Encoder
        self.enc1  = self._conv_block(in_channels, base_filters)
        self.pool1 = nn.MaxPool2d(2)

        self.enc2  = self._conv_block(base_filters, base_filters*2)
        self.pool2 = nn.MaxPool2d(2)

        self.enc3  = self._conv_block(base_filters*2, base_filters*4)
        self.pool3 = nn.MaxPool2d(2)

        self.enc4  = self._conv_block(base_filters*4, base_filters*8)
        self.pool4 = nn.MaxPool2d(2)

        # Bottleneck
        self.bottleneck = nn.Sequential(
            nn.Conv2d(base_filters*8, base_filters*16, kernel_size=3, padding=2, dilation=2),
            nn.BatchNorm2d(base_filters*16),
            nn.ReLU(),
            nn.Conv2d(base_filters*16, base_filters*16, kernel_size=3, padding=2, dilation=2),
            nn.BatchNorm2d(base_filters*16),
            nn.ReLU(),
        )

        # Decoder + skip connections
        self.up4  = nn.ConvTranspose2d(base_filters*16, base_filters*8,  kernel_size=2, stride=2)
        self.dec4 = self._conv_block(base_filters*8*2, base_filters*8)

        self.up3  = nn.ConvTranspose2d(base_filters*8,  base_filters*4,  kernel_size=2, stride=2)
        self.dec3 = self._conv_block(base_filters*4*2, base_filters*4)

        self.up2  = nn.ConvTranspose2d(base_filters*4,  base_filters*2,  kernel_size=2, stride=2)
        self.dec2 = self._conv_block(base_filters*2*2, base_filters*2)

        self.up1  = nn.ConvTranspose2d(base_filters*2,  base_filters,    kernel_size=2, stride=2)
        self.dec1 = self._conv_block(base_filters*2,   base_filters)

        # Output head
        self.conv_out    = nn.Conv2d(base_filters, 1, kernel_size=1)
        self.softplus_mu = nn.Softplus(beta=1.0, threshold=20.0)

        # Learnable log-dispersion parameter
        self.log_theta = nn.Parameter(torch.log(torch.tensor(5.0)))

        self.lr = lr
        self.weight_decay = weight_decay

    def _conv_block(self, in_ch, out_ch):
        return nn.Sequential(
            nn.Conv2d(in_ch, out_ch, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(),
            nn.Conv2d(out_ch, out_ch, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(),
        )

    def forward(self, x):
        # Encoder
        e1 = self.enc1(x);   p1 = self.pool1(e1)
        e2 = self.enc2(p1);  p2 = self.pool2(e2)
        e3 = self.enc3(p2);  p3 = self.pool3(e3)
        e4 = self.enc4(p3);  p4 = self.pool4(e4)

        # Bottleneck
        bn = self.bottleneck(p4)

        # Decoder
        d4 = self.up4(bn)
        d4 = self.dec4(torch.cat([d4, e4], dim=1))

        d3 = self.up3(d4)
        d3 = self.dec3(torch.cat([d3, e3], dim=1))

        d2 = self.up2(d3)
        d2 = self.dec2(torch.cat([d2, e2], dim=1))

        d1 = self.up1(d2)
        d1 = self.dec1(torch.cat([d1, e1], dim=1))

        # Output μ-map & θ
        mu_map = self.softplus_mu(self.conv_out(d1))
        mu     = mu_map.mean(dim=2).squeeze(1)

        # Stable θ: clamp log_theta then exp
        logt  = self.log_theta.clamp(min=math.log(1e-3), max=math.log(1e3))
        theta = logt.exp()

        return mu, theta

    def training_step(self, batch, batch_idx):
        if batch is None:
            return None
        x, y    = batch
        mu, θ   = self(x)
        mu      = mu.clamp(min=1e-3)
        p       = (θ / (θ + mu)).clamp(min=1e-6, max=1-1e-6)

        # manual NB log-prob
        lp = (
            torch.lgamma(y + θ)
          - torch.lgamma(θ)
          - torch.lgamma(y + 1)
          + θ * torch.log1p(-p)
          + y * torch.log(p)
        )
        loss = -lp.mean()

        # --- original on-metal distribution (commented) ---
        # from torch.distributions import NegativeBinomial
        # nb = NegativeBinomial(total_count=θ, probs=p)
        # loss = -nb.log_prob(y).mean()

        self.log("train_loss", loss)
        return loss

    def validation_step(self, batch, batch_idx):
        if batch is None:
            return None        
        x, y    = batch
        mu, θ   = self(x)
        mu      = mu.clamp(min=1e-3)
        p       = (θ / (θ + mu)).clamp(min=1e-6, max=1-1e-6)

        # metrics
        self.log("val_mu_min",    mu.min(),    prog_bar=True)
        self.log("val_mu_max",    mu.max(),    prog_bar=True)
        self.log("val_probs_min", p.min(),     prog_bar=True)
        self.log("val_probs_max", p.max(),     prog_bar=True)

        # manual NB log-prob on CPU for extra safety
        y_cpu    = y.detach().cpu()
        μ_cpu    = mu.detach().cpu()
        θ_cpu    = θ.detach().cpu()
        p_cpu    = p.detach().cpu()
        lp_cpu   = (
            torch.lgamma(y_cpu + θ_cpu)
          - torch.lgamma(θ_cpu)
          - torch.lgamma(y_cpu + 1)
          + θ_cpu * torch.log1p(-p_cpu)
          + y_cpu * torch.log(p_cpu)
        )
        loss_cpu = -lp_cpu.mean()
        loss     = loss_cpu.to(mu.device)

        # --- original on-metal distribution (commented) ---
        # nb_cpu = NegativeBinomial(total_count=θ_cpu, probs=p_cpu)
        # loss   = -(nb_cpu.log_prob(y_cpu).mean()).to(mu.device)

        self.log("val_loss", loss, prog_bar=True)
        return loss

    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        if batch is None:
            return None
        x, _       = batch
        mu, theta  = self(x)
        return {"mu": mu, "theta": theta}

    def on_after_backward(self):
        # clamp log_theta after grad update
        self.log_theta.data.clamp_(min=math.log(1e-3), max=math.log(1e3))

    def configure_optimizers(self):
        return torch.optim.AdamW(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)

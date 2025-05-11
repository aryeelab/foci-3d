import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import pytorch_lightning as pl

class FootprintDataset(Dataset):
    def __init__(self, regions, accessor_fn):
        """
        regions: list of (chrom, start, end) tuples
        accessor_fn: function(chrom, start, end) -> (
          footprint: np.ndarray of shape (S, W),
          total_signal: ignored,
          procap: np.ndarray of shape (W,)
        )
        The footprint is automatically zero-padded to the nearest multiple of 16 in height.
        """
        self.regions   = regions
        self.accessor  = accessor_fn

    def __len__(self):
        return len(self.regions)

    def __getitem__(self, idx):
        chrom, start, end = self.regions[idx]
        fp, _, procap   = self.accessor(chrom, start, end)

        # convert pandas to numpy if needed
        if hasattr(fp, "to_numpy"):
            fp = fp.to_numpy()
        if hasattr(procap, "to_numpy"):
            procap = procap.to_numpy()

        # ensure height is divisible by 16 for U-Net pooling
        S, W = fp.shape
        target_H = ((S + 15) // 16) * 16
        pad_total = target_H - S
        if pad_total > 0:
            pad_top = pad_total // 2
            pad_bottom = pad_total - pad_top
            fp = np.pad(fp, ((pad_top, pad_bottom), (0, 0)), mode='constant')

        # Log transform for training stability
        fp = np.sign(fp) * np.log1p(np.abs(fp))

        # input tensor: single-channel image (1, H, W)
        x = torch.from_numpy(fp).float().unsqueeze(0)

        # target tensor: per-position signal (W,)
        # NegativeBinomial expects integer counts, so we round and cast
        y = torch.from_numpy(procap).float().round().long()

        return x, y

class FootprintDataModule(pl.LightningDataModule):
    def __init__(self, train_regions, val_regions, test_regions, accessor_fn,
                 batch_size=16, num_workers=4):
        super().__init__()
        self.train_regions   = train_regions
        self.val_regions     = val_regions
        self.test_regions    = test_regions
        self.accessor_fn     = accessor_fn
        self.batch_size      = batch_size
        self.num_workers     = num_workers

    def setup(self, stage=None):
        self.train_ds = FootprintDataset(self.train_regions, self.accessor_fn)
        self.val_ds   = FootprintDataset(self.val_regions,   self.accessor_fn)
        self.test_ds  = FootprintDataset(self.test_regions,  self.accessor_fn)  
        
    def train_dataloader(self):
        return DataLoader(
            self.train_ds,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
            persistent_workers=True
        )

    def val_dataloader(self):
        return DataLoader(
            self.val_ds,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            persistent_workers=True
        )

    def test_dataloader(self):
        return DataLoader(
            self.test_ds,
            batch_size=1,
            shuffle=False,
            num_workers=self.num_workers,
            persistent_workers=True
        )
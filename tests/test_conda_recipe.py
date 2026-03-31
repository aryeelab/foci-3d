import unittest
from pathlib import Path


class CondaRecipeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.repo_root = Path(__file__).resolve().parent.parent
        self.recipe_text = (self.repo_root / "conda-recipe" / "meta.yaml").read_text()

    def test_recipe_has_release_sha256(self) -> None:
        self.assertNotIn("REPLACE_WITH_RELEASE_TARBALL_SHA256", self.recipe_text)
        self.assertIn("sha256:", self.recipe_text)

    def test_recipe_has_run_exports_pin(self) -> None:
        self.assertIn("run_exports:", self.recipe_text)
        self.assertIn('{{ pin_subpackage(name, max_pin="x.x") }}', self.recipe_text)

    def test_recipe_points_at_github_release_tarball(self) -> None:
        self.assertIn(
            "https://github.com/aryeelab/foci-3d/archive/refs/tags/v{{ version }}.tar.gz",
            self.recipe_text,
        )


if __name__ == "__main__":
    unittest.main()

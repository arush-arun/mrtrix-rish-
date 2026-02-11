"""
Unit tests for voxel-wise GLM module.

Tests the core statistical computations following MRtrix3 design patterns.
"""

import numpy as np
import pytest

from src.qc.glm import (
    Hypothesis,
    Partition,
    TestFixedHomoscedastic,
    TestFixedHeteroscedastic,
    check_design,
    create_design_matrix,
    create_site_contrast,
)


class TestHypothesis:
    """Tests for Hypothesis class."""

    def test_t_test_contrast(self):
        """Test single-row contrast (t-test)."""
        contrast = np.array([1, -1, 0, 0])
        h = Hypothesis(contrast, index=0, name="site_diff")

        assert h.cols == 4
        assert h.rank == 1
        assert h.is_F is False
        assert h.name == "site_diff"
        assert h.index == 0

    def test_f_test_contrast(self):
        """Test multi-row contrast (F-test)."""
        contrast = np.array([
            [1, 0, 0, 0],
            [0, 1, 0, 0]
        ])
        h = Hypothesis(contrast, index=1)

        assert h.cols == 4
        assert h.rank == 2
        assert h.is_F is True
        assert h.name == "F1"

    def test_default_naming(self):
        """Test default hypothesis naming."""
        h_t = Hypothesis(np.array([1, -1]))
        h_f = Hypothesis(np.array([[1, 0], [0, 1]]), index=2)

        assert h_t.name == "t0"
        assert h_f.name == "F2"

    def test_partition_simple(self):
        """Test design matrix partitioning."""
        # Design: [intercept, site, covariate]
        n = 10
        design = np.column_stack([
            np.ones(n),
            np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1]),  # site
            np.random.randn(n)  # covariate
        ])

        # Contrast for site effect
        contrast = np.array([0, 1, 0])  # Test site coefficient
        h = Hypothesis(contrast)

        part = h.partition(design)

        assert isinstance(part, Partition)
        assert part.X.shape[1] == 1  # Site column
        assert part.Z.shape[1] == 2  # Intercept + covariate
        assert part.Rz.shape == (n, n)
        assert part.Hz.shape == (n, n)

    def test_partition_caching(self):
        """Test that partition is cached."""
        design = np.random.randn(10, 3)
        contrast = np.array([1, 0, 0])
        h = Hypothesis(contrast)

        part1 = h.partition(design)
        part2 = h.partition(design)

        assert part1 is part2  # Same object (cached)


class TestDesignMatrix:
    """Tests for design matrix creation."""

    def test_simple_two_sites(self):
        """Test design matrix with two sites."""
        site_labels = ["A", "A", "A", "B", "B", "B"]
        design, names = create_design_matrix(site_labels)

        assert design.shape == (6, 2)  # intercept + 1 site dummy
        assert names == ["intercept", "site_B"]
        assert np.allclose(design[:, 0], 1.0)  # intercept
        assert np.allclose(design[3:, 1], 1.0)  # site B
        assert np.allclose(design[:3, 1], 0.0)  # site A (reference)

    def test_three_sites(self):
        """Test design matrix with three sites."""
        site_labels = ["X", "X", "Y", "Y", "Z", "Z"]
        design, names = create_design_matrix(site_labels)

        assert design.shape == (6, 3)  # intercept + 2 dummies
        assert names == ["intercept", "site_Y", "site_Z"]

    def test_with_covariates(self):
        """Test design matrix with covariates."""
        site_labels = ["A", "A", "B", "B"]
        covariates = {
            "age": [25, 30, 35, 40],
            "sex": [0, 1, 0, 1]
        }
        design, names = create_design_matrix(site_labels, covariates)

        assert design.shape == (4, 4)  # intercept + site + 2 covariates
        assert names == ["intercept", "site_B", "age", "sex"]

    def test_covariate_standardization(self):
        """Test that covariates are standardized by default."""
        site_labels = ["A", "A", "B", "B"]
        covariates = {"age": [20, 30, 40, 50]}
        design, _ = create_design_matrix(site_labels, covariates, standardize_covariates=True)

        # Age column should be standardized (mean ~0, std ~1)
        age_col = design[:, -1]
        assert np.abs(age_col.mean()) < 1e-10
        assert np.abs(age_col.std() - 1.0) < 0.1


class TestCheckDesign:
    """Tests for design matrix validation."""

    def test_full_rank(self):
        """Test full rank design matrix."""
        design = np.array([
            [1, 0, 0.5],
            [1, 0, 0.8],
            [1, 1, 0.3],
            [1, 1, 0.7]
        ])
        rank, cond = check_design(design)

        assert rank == 3
        assert cond < 100

    def test_rank_deficient(self):
        """Test rank-deficient design matrix."""
        # Duplicate column
        design = np.array([
            [1, 1, 0],
            [1, 1, 0],
            [1, 1, 1],
            [1, 1, 1]
        ])
        rank, cond = check_design(design)

        assert rank < 3  # Not full rank

    def test_high_condition(self):
        """Test high condition number warning case."""
        # Highly correlated columns
        design = np.array([
            [1, 1.0001, 0],
            [1, 1.0002, 0],
            [1, 1.0003, 1],
            [1, 1.0004, 1]
        ])
        rank, cond = check_design(design)

        assert cond > 100  # High condition number


class TestSiteContrast:
    """Tests for site effect contrast creation."""

    def test_two_sites(self):
        """Test contrast for two sites."""
        h = create_site_contrast(n_sites=2)

        assert h.is_F is False  # Single row = t-test
        assert h.cols == 2  # intercept + 1 site dummy
        assert h.name == "site_effect"

    def test_three_sites(self):
        """Test contrast for three sites."""
        h = create_site_contrast(n_sites=3)

        assert h.is_F is True  # Multiple rows = F-test
        assert h.cols == 3  # intercept + 2 site dummies
        assert h.rank == 2


class TestTestFixedHomoscedastic:
    """Tests for homoscedastic GLM test."""

    @pytest.fixture
    def synthetic_data(self):
        """Create synthetic data with known site effect."""
        np.random.seed(42)
        n_per_site = 20
        n_voxels = 100

        # Two sites with different means
        site_A = np.random.randn(n_per_site, n_voxels) + 0.0
        site_B = np.random.randn(n_per_site, n_voxels) + 1.0  # Effect size = 1

        data = np.vstack([site_A, site_B])
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site

        return data, site_labels

    def test_solve_betas(self, synthetic_data):
        """Test beta coefficient estimation."""
        data, site_labels = synthetic_data
        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)

        test = TestFixedHomoscedastic(data, design, [hypothesis])
        betas = test.solve_betas()

        assert betas.shape == (2, 100)  # 2 predictors × 100 voxels

        # Site B effect should be approximately 1.0
        site_effect = betas[1, :].mean()
        assert 0.5 < site_effect < 1.5

    def test_residuals(self, synthetic_data):
        """Test residual computation."""
        data, site_labels = synthetic_data
        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)

        test = TestFixedHomoscedastic(data, design, [hypothesis])
        residuals = test.residuals

        assert residuals.shape == data.shape
        # Residuals should have mean ~0 per column
        assert np.abs(residuals.mean(axis=0)).max() < 0.1

    def test_f_statistic(self, synthetic_data):
        """Test F-statistic computation."""
        data, site_labels = synthetic_data
        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)

        test = TestFixedHomoscedastic(data, design, [hypothesis])
        outputs = test()

        assert len(outputs) == 1
        f_stat = outputs[0].statistic

        assert f_stat.shape == (100,)
        # With effect size ~1 and n=40, F should be notably larger than 1
        # (F ~ n * d^2 / 2 where d=1, n=40, so expect F ~ 20, but variance reduces this)
        assert f_stat.mean() > 3.0  # Reasonable threshold for effect size 1

    def test_effect_size(self, synthetic_data):
        """Test effect size computation."""
        data, site_labels = synthetic_data
        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)

        test = TestFixedHomoscedastic(data, design, [hypothesis])
        outputs = test()

        # For t-test (2 sites), we get effect sizes
        effect = outputs[0].effect_size
        assert effect is not None
        # Effect should be ~1.0
        assert 0.5 < effect.mean() < 1.5

    def test_permutation(self, synthetic_data):
        """Test permutation inference."""
        data, site_labels = synthetic_data
        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)

        test = TestFixedHomoscedastic(data, design, [hypothesis])

        # Observed
        observed = test()[0].statistic

        # Permuted
        perm_indices = np.random.permutation(len(site_labels))
        permuted = test(perm_indices)[0].statistic

        # Permuted F should be smaller (destroys signal)
        assert permuted.mean() < observed.mean()

    def test_no_site_effect(self):
        """Test with no site effect (null case)."""
        np.random.seed(42)
        n_per_site = 20
        n_voxels = 50

        # Both sites have same distribution
        data = np.random.randn(2 * n_per_site, n_voxels)
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site

        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)

        test = TestFixedHomoscedastic(data, design, [hypothesis])
        outputs = test()

        f_stat = outputs[0].statistic
        # F should be close to 1 (null distribution mean)
        assert f_stat.mean() < 5.0


class TestTestFixedHeteroscedastic:
    """Tests for heteroscedastic GLM test."""

    def test_variance_groups(self):
        """Test heteroscedastic test with different variance groups."""
        np.random.seed(42)
        n_per_site = 20
        n_voxels = 50

        # Site A: low variance, Site B: high variance
        site_A = np.random.randn(n_per_site, n_voxels) * 0.5 + 0.0
        site_B = np.random.randn(n_per_site, n_voxels) * 2.0 + 1.0

        data = np.vstack([site_A, site_B])
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site
        variance_groups = np.array([0] * n_per_site + [1] * n_per_site)

        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)

        test = TestFixedHeteroscedastic(data, design, [hypothesis], variance_groups)
        outputs = test()

        assert len(outputs) == 1
        # Should produce valid statistics
        assert not np.any(np.isnan(outputs[0].statistic))


class TestGLMIntegration:
    """Integration tests for complete GLM workflow."""

    def test_three_sites_with_covariates(self):
        """Test three-site analysis with covariates."""
        np.random.seed(42)
        n_per_site = 15
        n_voxels = 30

        # Three sites with different effects
        site_A = np.random.randn(n_per_site, n_voxels)
        site_B = np.random.randn(n_per_site, n_voxels) + 0.5
        site_C = np.random.randn(n_per_site, n_voxels) + 1.0

        data = np.vstack([site_A, site_B, site_C])
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site + ["C"] * n_per_site

        # Add age covariate
        ages = np.concatenate([
            np.random.uniform(20, 40, n_per_site),
            np.random.uniform(30, 50, n_per_site),
            np.random.uniform(40, 60, n_per_site)
        ])
        covariates = {"age": ages.tolist()}

        design, names = create_design_matrix(site_labels, covariates)
        hypothesis = create_site_contrast(n_sites=3, n_covariates=1)

        test = TestFixedHomoscedastic(data, design, [hypothesis])
        outputs = test()

        # Should detect site effect
        assert outputs[0].statistic.mean() > 1.0
        assert hypothesis.is_F  # F-test for 3 sites


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

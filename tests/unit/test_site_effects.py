"""
Unit tests for site effect statistical testing module.

Tests permutation framework, FDR correction, and effect size computation.
"""

import numpy as np
import pytest

from src.qc.site_effects import (
    Shuffle,
    Shuffler,
    fdr_correction,
    permutation_p_values,
    compute_partial_eta_squared,
    compute_cohens_f,
)


class TestShuffle:
    """Tests for Shuffle dataclass."""

    def test_shuffle_creation(self):
        """Test creating a shuffle object."""
        indices = np.array([2, 0, 1, 3])
        s = Shuffle(index=1, data=indices)

        assert s.index == 1
        assert np.array_equal(s.data, indices)


class TestShuffler:
    """Tests for Shuffler permutation generator."""

    def test_basic_permutations(self):
        """Test basic permutation generation."""
        shuffler = Shuffler(n_subjects=10, n_permutations=100, seed=42)

        assert len(shuffler) == 100  # Including identity

    def test_identity_first(self):
        """Test that first permutation is identity."""
        shuffler = Shuffler(n_subjects=5, n_permutations=10, seed=42)
        shuffles = list(shuffler)

        identity = shuffles[0]
        assert identity.index == 0
        assert np.array_equal(identity.data, np.arange(5))

    def test_no_duplicates(self):
        """Test that permutations are unique."""
        shuffler = Shuffler(n_subjects=10, n_permutations=50, seed=42)
        shuffles = list(shuffler)

        # Convert to tuples for comparison
        perm_tuples = [tuple(s.data) for s in shuffles]
        assert len(perm_tuples) == len(set(perm_tuples))

    def test_permutation_validity(self):
        """Test that each permutation is a valid reordering."""
        shuffler = Shuffler(n_subjects=8, n_permutations=20, seed=42)

        for shuffle in shuffler:
            # Should contain all indices 0 to n-1
            assert set(shuffle.data) == set(range(8))
            assert len(shuffle.data) == 8

    def test_seed_reproducibility(self):
        """Test that seed produces reproducible results."""
        shuffler1 = Shuffler(n_subjects=10, n_permutations=20, seed=123)
        shuffler2 = Shuffler(n_subjects=10, n_permutations=20, seed=123)

        perms1 = [tuple(s.data) for s in shuffler1]
        perms2 = [tuple(s.data) for s in shuffler2]

        assert perms1 == perms2

    def test_different_seeds(self):
        """Test that different seeds produce different results."""
        shuffler1 = Shuffler(n_subjects=10, n_permutations=20, seed=123)
        shuffler2 = Shuffler(n_subjects=10, n_permutations=20, seed=456)

        perms1 = [tuple(s.data) for s in shuffler1]
        perms2 = [tuple(s.data) for s in shuffler2]

        # At least some should differ (except identity)
        assert perms1[1:] != perms2[1:]

    def test_iterator_reset(self):
        """Test iterator reset functionality."""
        shuffler = Shuffler(n_subjects=5, n_permutations=10, seed=42)

        first_pass = list(shuffler)
        shuffler.reset()
        second_pass = list(shuffler)

        # Should produce same sequence
        for s1, s2 in zip(first_pass, second_pass):
            assert np.array_equal(s1.data, s2.data)

    def test_exchangeability_blocks(self):
        """Test permutation with exchangeability blocks."""
        # Block structure: subjects 0-2 in block 0, subjects 3-5 in block 1
        eb = np.array([0, 0, 0, 1, 1, 1])
        shuffler = Shuffler(
            n_subjects=6,
            n_permutations=50,
            exchangeability_blocks=eb,
            seed=42
        )

        for shuffle in shuffler:
            if shuffle.index == 0:
                continue  # Skip identity

            # Check that permutation respects blocks
            perm = shuffle.data
            block0_positions = set(perm[:3])
            block1_positions = set(perm[3:])

            # Each block should only contain its own indices
            assert block0_positions == {0, 1, 2}
            assert block1_positions == {3, 4, 5}


class TestFDRCorrection:
    """Tests for FDR correction."""

    def test_all_significant(self):
        """Test when all p-values are very small."""
        p_values = np.array([0.001, 0.002, 0.003, 0.004, 0.005])
        q_values, threshold, significant = fdr_correction(p_values, alpha=0.05)

        assert np.all(significant)
        assert threshold > 0

    def test_none_significant(self):
        """Test when all p-values are large."""
        p_values = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        q_values, threshold, significant = fdr_correction(p_values, alpha=0.05)

        assert not np.any(significant)
        assert threshold == 0.0

    def test_mixed_significance(self):
        """Test with mix of significant and non-significant."""
        p_values = np.array([0.001, 0.01, 0.05, 0.1, 0.5])
        q_values, threshold, significant = fdr_correction(p_values, alpha=0.05)

        # Smallest p-values should be significant
        assert significant[0]
        # Largest should not be
        assert not significant[-1]

    def test_q_values_ordering(self):
        """Test that q-values preserve ordering."""
        p_values = np.array([0.01, 0.05, 0.001, 0.1, 0.03])
        q_values, _, _ = fdr_correction(p_values, alpha=0.05)

        # Smaller p should have smaller q
        p_order = np.argsort(p_values)
        q_order = np.argsort(q_values)
        assert np.array_equal(p_order, q_order)

    def test_with_nan(self):
        """Test handling of NaN values."""
        p_values = np.array([0.01, np.nan, 0.05, 0.1])
        q_values, threshold, significant = fdr_correction(p_values, alpha=0.05)

        # NaN should remain NaN
        assert np.isnan(q_values[1])
        # Other values should be processed
        assert not np.isnan(q_values[0])

    def test_benjamini_yekutieli(self):
        """Test BY correction (more conservative)."""
        p_values = np.array([0.001, 0.01, 0.02, 0.03, 0.04])

        q_bh, thresh_bh, sig_bh = fdr_correction(p_values, alpha=0.05, method="bh")
        q_by, thresh_by, sig_by = fdr_correction(p_values, alpha=0.05, method="by")

        # BY should be more conservative (fewer significant)
        assert sig_by.sum() <= sig_bh.sum()


class TestPermutationPValues:
    """Tests for permutation-based p-value computation."""

    def test_extreme_observed(self):
        """Test with extreme observed value."""
        observed = np.array([100.0])
        null = np.random.randn(1000, 1) * 10  # Max ~30

        p_values = permutation_p_values(observed, null, tail="right")

        # Should be very small
        assert p_values[0] < 0.01

    def test_null_observed(self):
        """Test with observed in null range."""
        observed = np.array([0.0])
        null = np.random.randn(1000, 1)

        p_values = permutation_p_values(observed, null, tail="right")

        # Should be around 0.5
        assert 0.3 < p_values[0] < 0.7

    def test_multiple_elements(self):
        """Test with multiple elements."""
        observed = np.array([10.0, 0.0, -5.0])
        null = np.random.randn(1000, 3) * 2

        p_values = permutation_p_values(observed, null, tail="right")

        assert len(p_values) == 3
        assert p_values[0] < p_values[1]  # 10 more extreme than 0

    def test_two_tailed(self):
        """Test two-tailed p-values."""
        observed = np.array([-10.0, 10.0])
        null = np.random.randn(1000, 2) * 2

        p_right = permutation_p_values(observed, null, tail="right")
        p_two = permutation_p_values(observed, null, tail="two")

        # Both extreme values should have small two-tailed p
        assert p_two[0] < 0.01
        assert p_two[1] < 0.01

    def test_minimum_p_value(self):
        """Test that p-value is at least 1/(n_perms+1)."""
        observed = np.array([1000.0])  # Extremely large
        null = np.zeros((100, 1))

        p_values = permutation_p_values(observed, null, tail="right")

        # Minimum is 1/101
        assert p_values[0] == pytest.approx(1 / 101, rel=1e-10)


class TestEffectSize:
    """Tests for effect size computations."""

    def test_partial_eta_squared_strong_effect(self):
        """Test eta-squared with strong site effect."""
        np.random.seed(42)

        # Two sites with very different means
        n_per_site = 50
        n_voxels = 10

        site_A = np.random.randn(n_per_site, n_voxels)
        site_B = np.random.randn(n_per_site, n_voxels) + 5.0  # Large effect

        data = np.vstack([site_A, site_B])
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site

        eta_sq = compute_partial_eta_squared(data, site_labels)

        # Should be large (>0.5)
        assert eta_sq.mean() > 0.5

    def test_partial_eta_squared_no_effect(self):
        """Test eta-squared with no site effect."""
        np.random.seed(42)

        n_per_site = 50
        n_voxels = 10

        # Same distribution for both sites
        data = np.random.randn(2 * n_per_site, n_voxels)
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site

        eta_sq = compute_partial_eta_squared(data, site_labels)

        # Should be small (<0.1)
        assert eta_sq.mean() < 0.1

    def test_partial_eta_squared_range(self):
        """Test that eta-squared is in [0, 1]."""
        np.random.seed(42)

        data = np.random.randn(40, 20)
        site_labels = ["A"] * 20 + ["B"] * 20

        eta_sq = compute_partial_eta_squared(data, site_labels)

        assert np.all(eta_sq >= 0)
        assert np.all(eta_sq <= 1)

    def test_cohens_f_conversion(self):
        """Test conversion from eta-squared to Cohen's f."""
        eta_sq = np.array([0.01, 0.06, 0.14, 0.25])

        f = compute_cohens_f(eta_sq)

        # Small, medium, large effects (Cohen's conventions)
        # f = 0.10 is small, 0.25 is medium, 0.40 is large
        assert f[0] <= 0.11  # Small (~0.1005)
        assert 0.20 < f[1] < 0.30  # Small-medium (~0.253)
        assert 0.35 < f[2] < 0.45  # Medium-large (~0.404)
        assert f[3] > 0.5  # Large (~0.577)

    def test_cohens_f_formula(self):
        """Test Cohen's f formula: f = sqrt(eta / (1 - eta))."""
        eta_sq = np.array([0.04, 0.09, 0.16])

        f = compute_cohens_f(eta_sq)

        expected = np.sqrt(eta_sq / (1 - eta_sq))
        assert np.allclose(f, expected)


class TestSiteEffectIntegration:
    """Integration tests for site effect testing workflow."""

    def test_synthetic_site_effect(self):
        """Test complete workflow with synthetic data."""
        from src.qc.glm import (
            TestFixedHomoscedastic,
            create_design_matrix,
            create_site_contrast,
        )

        np.random.seed(42)
        n_per_site = 30
        n_voxels = 50

        # Strong site effect
        site_A = np.random.randn(n_per_site, n_voxels)
        site_B = np.random.randn(n_per_site, n_voxels) + 1.5

        data = np.vstack([site_A, site_B])
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site

        # Create model
        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)
        test = TestFixedHomoscedastic(data, design, [hypothesis])

        # Observed statistics
        observed = test()[0].statistic

        # Permutation null
        n_perms = 100
        null_dist = []
        for i in range(n_perms):
            perm = np.random.permutation(len(site_labels))
            null_dist.append(test(perm)[0].statistic)
        null_dist = np.array(null_dist)

        # P-values
        p_perm = permutation_p_values(observed, null_dist, tail="right")

        # FDR correction
        _, _, significant = fdr_correction(p_perm, alpha=0.05)

        # Should detect many significant voxels
        assert significant.mean() > 0.5

    def test_no_site_effect_null(self):
        """Test that null case produces uniform p-values."""
        from src.qc.glm import (
            TestFixedHomoscedastic,
            create_design_matrix,
            create_site_contrast,
        )

        np.random.seed(42)
        n_per_site = 30
        n_voxels = 100

        # No site effect
        data = np.random.randn(2 * n_per_site, n_voxels)
        site_labels = ["A"] * n_per_site + ["B"] * n_per_site

        design, _ = create_design_matrix(site_labels)
        hypothesis = create_site_contrast(n_sites=2)
        test = TestFixedHomoscedastic(data, design, [hypothesis])

        observed = test()[0].statistic

        # Permutation null
        n_perms = 200
        null_dist = []
        for i in range(n_perms):
            perm = np.random.permutation(len(site_labels))
            null_dist.append(test(perm)[0].statistic)
        null_dist = np.array(null_dist)

        p_perm = permutation_p_values(observed, null_dist, tail="right")

        # FDR correction
        _, _, significant = fdr_correction(p_perm, alpha=0.05)

        # Should have few significant voxels (around 5%)
        assert significant.mean() < 0.15


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

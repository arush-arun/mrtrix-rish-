"""
Unit tests for RISH-GLM harmonization module.

Tests the joint GLM-based estimation of site effects and covariates.
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from src.core.rish_glm import (
    RISHGLMResult,
    build_rish_glm_design,
    save_rish_glm_model,
    load_rish_glm_model,
)


class TestBuildDesignMatrix:
    """Tests for design matrix construction."""

    def test_two_sites_no_covariates(self):
        """Two sites, no covariates -> 2-column indicator matrix."""
        labels = ["A", "A", "A", "B", "B"]
        design, names, site_map, means, stds = build_rish_glm_design(labels)

        assert design.shape == (5, 2)
        assert names == ["site_A", "site_B"]
        assert site_map == {"A": 0, "B": 1}
        assert means == {}
        assert stds == {}

        # Check indicators
        expected = np.array([
            [1, 0],
            [1, 0],
            [1, 0],
            [0, 1],
            [0, 1],
        ], dtype=np.float64)
        np.testing.assert_array_equal(design, expected)

    def test_three_sites_no_covariates(self):
        """Three sites -> 3-column indicator matrix."""
        labels = ["X", "Y", "X", "Z", "Y", "Z"]
        design, names, site_map, _, _ = build_rish_glm_design(labels)

        assert design.shape == (6, 3)
        assert names == ["site_X", "site_Y", "site_Z"]
        assert site_map == {"X": 0, "Y": 1, "Z": 2}

        # Each row should have exactly one 1
        assert np.all(design.sum(axis=1) == 1.0)

    def test_no_intercept(self):
        """Design matrix should NOT have an intercept column."""
        labels = ["A", "A", "B", "B"]
        design, names, _, _, _ = build_rish_glm_design(labels)

        # No column of all ones
        for j in range(design.shape[1]):
            assert not np.all(design[:, j] == 1.0), \
                f"Column {j} ({names[j]}) is an intercept"

    def test_with_covariates(self):
        """Two sites + 2 covariates -> 4-column matrix."""
        labels = ["A", "A", "B", "B"]
        covariates = {
            "age": [30.0, 40.0, 35.0, 45.0],
            "sex": [0.0, 1.0, 0.0, 1.0],
        }
        design, names, site_map, means, stds = build_rish_glm_design(
            labels, covariates
        )

        assert design.shape == (4, 4)
        assert names == ["site_A", "site_B", "age", "sex"]
        assert site_map == {"A": 0, "B": 1}

        # Covariates should be z-scored
        age_col = design[:, 2]
        assert abs(age_col.mean()) < 1e-10, "Age not zero-mean"
        assert abs(age_col.std() - 1.0) < 0.1, "Age not unit-variance"

        assert "age" in means
        assert "sex" in means
        assert "age" in stds
        assert "sex" in stds

    def test_covariates_sorted(self):
        """Covariates should be in sorted order."""
        labels = ["A", "B"]
        covariates = {"zebra": [1.0, 2.0], "alpha": [3.0, 4.0]}
        _, names, _, _, _ = build_rish_glm_design(labels, covariates)

        assert names == ["site_A", "site_B", "alpha", "zebra"]

    def test_covariate_length_mismatch(self):
        """Covariate with wrong length should raise ValueError."""
        labels = ["A", "B", "A"]
        covariates = {"age": [30.0, 40.0]}  # 2 values, 3 subjects

        with pytest.raises(ValueError, match="expected 3"):
            build_rish_glm_design(labels, covariates)

    def test_constant_covariate(self):
        """Constant covariate should not cause division by zero."""
        labels = ["A", "B", "A"]
        covariates = {"group": [1.0, 1.0, 1.0]}
        design, _, _, _, stds = build_rish_glm_design(labels, covariates)

        # std defaults to 1.0 for constant covariates
        assert stds["group"] == 1.0
        assert np.all(np.isfinite(design))


class TestGLMFittingSynthetic:
    """Test GLM fitting with synthetic data (no MRtrix3 dependency)."""

    def test_beta_recovery_no_covariates(self):
        """Verify betas recover true site means from synthetic data."""
        np.random.seed(42)
        n_ref, n_tar = 10, 8
        n_voxels = 50
        true_ref_mean = 1.0
        true_tar_mean = 0.8

        # Simulate RISH data
        ref_data = np.random.normal(true_ref_mean, 0.05, (n_ref, n_voxels))
        tar_data = np.random.normal(true_tar_mean, 0.05, (n_tar, n_voxels))
        data = np.vstack([ref_data, tar_data])

        labels = ["ref"] * n_ref + ["tar"] * n_tar
        design, _, site_map, _, _ = build_rish_glm_design(labels)

        # Fit GLM
        betas, _, _, _ = np.linalg.lstsq(design, data, rcond=None)

        # Check betas recover true means
        ref_col = site_map["ref"]
        tar_col = site_map["tar"]
        np.testing.assert_allclose(
            betas[ref_col].mean(), true_ref_mean, atol=0.05
        )
        np.testing.assert_allclose(
            betas[tar_col].mean(), true_tar_mean, atol=0.05
        )

    def test_scale_factor_from_betas(self):
        """Scale factor β_ref / β_tar should recover true ratio."""
        np.random.seed(42)
        n_ref, n_tar = 20, 20
        n_voxels = 100
        true_ref_mean = 1.0
        true_tar_mean = 0.5  # Expected scale = 2.0

        ref_data = np.random.normal(true_ref_mean, 0.02, (n_ref, n_voxels))
        tar_data = np.random.normal(true_tar_mean, 0.02, (n_tar, n_voxels))
        data = np.vstack([ref_data, tar_data])

        labels = ["ref"] * n_ref + ["tar"] * n_tar
        design, _, site_map, _, _ = build_rish_glm_design(labels)

        betas, _, _, _ = np.linalg.lstsq(design, data, rcond=None)

        # Scale factor
        scale = betas[site_map["ref"]] / np.maximum(
            betas[site_map["tar"]], 0.01
        )
        np.testing.assert_allclose(scale.mean(), 2.0, atol=0.1)

    def test_covariate_corrects_confound(self):
        """Including covariates should correctly separate site from age."""
        np.random.seed(42)
        n_ref, n_tar = 15, 15
        n_voxels = 80
        true_ref_mean = 1.0
        true_tar_mean = 0.8
        age_effect = 0.01  # per year

        # Ages: reference site younger, target site older (confounded)
        ages_ref = np.random.normal(30, 5, n_ref)
        ages_tar = np.random.normal(50, 5, n_tar)
        ages = np.concatenate([ages_ref, ages_tar])
        pop_mean_age = ages.mean()

        # Data: site mean + age effect + noise
        ref_data = true_ref_mean + age_effect * ages_ref[:, None] + \
            np.random.normal(0, 0.02, (n_ref, n_voxels))
        tar_data = true_tar_mean + age_effect * ages_tar[:, None] + \
            np.random.normal(0, 0.02, (n_tar, n_voxels))
        data = np.vstack([ref_data, tar_data])

        labels = ["ref"] * n_ref + ["tar"] * n_tar
        covariates = {"age": ages.tolist()}

        # Fit WITH covariates
        design_cov, _, site_map_cov, _, _ = build_rish_glm_design(
            labels, covariates
        )
        betas_cov, _, _, _ = np.linalg.lstsq(design_cov, data, rcond=None)

        # Fit WITHOUT covariates
        design_no, _, site_map_no, _, _ = build_rish_glm_design(labels)
        betas_no, _, _, _ = np.linalg.lstsq(design_no, data, rcond=None)

        # β_site represents site mean at z_age=0, i.e., at population mean age
        expected_ref = true_ref_mean + age_effect * pop_mean_age
        expected_tar = true_tar_mean + age_effect * pop_mean_age
        np.testing.assert_allclose(
            betas_cov[site_map_cov["ref"]].mean(), expected_ref, atol=0.05
        )
        np.testing.assert_allclose(
            betas_cov[site_map_cov["tar"]].mean(), expected_tar, atol=0.05
        )

        # Key test: scale factor WITH covariates correctly reflects
        # the pure site difference (not confounded by age)
        scale_cov = (betas_cov[site_map_cov["ref"]] /
                     np.maximum(betas_cov[site_map_cov["tar"]], 0.01))
        expected_scale = expected_ref / expected_tar
        np.testing.assert_allclose(
            scale_cov.mean(), expected_scale, atol=0.05
        )

        # WITHOUT covariates, scale factor is different because age
        # confound makes the two site means closer together
        scale_no = (betas_no[site_map_no["ref"]] /
                    np.maximum(betas_no[site_map_no["tar"]], 0.01))
        # ref mean ≈ 1.0 + 0.01*30 = 1.3, tar mean ≈ 0.8 + 0.01*50 = 1.3
        # so without covariates, scale ≈ 1.0 (age confound masks site diff)
        assert abs(scale_no.mean() - expected_scale) > 0.05, \
            "Without covariates, scale should differ from true value"

    def test_three_sites(self):
        """Three-site model should produce correct betas."""
        np.random.seed(42)
        n_a, n_b, n_c = 10, 10, 10
        n_voxels = 30
        means = {"A": 1.0, "B": 0.8, "C": 0.6}

        data = np.vstack([
            np.random.normal(means["A"], 0.02, (n_a, n_voxels)),
            np.random.normal(means["B"], 0.02, (n_b, n_voxels)),
            np.random.normal(means["C"], 0.02, (n_c, n_voxels)),
        ])
        labels = ["A"] * n_a + ["B"] * n_b + ["C"] * n_c

        design, _, site_map, _, _ = build_rish_glm_design(labels)
        betas, _, _, _ = np.linalg.lstsq(design, data, rcond=None)

        for site, true_mean in means.items():
            np.testing.assert_allclose(
                betas[site_map[site]].mean(), true_mean, atol=0.05
            )


class TestModelSaveLoad:
    """Test model serialization round-trip."""

    def test_save_load_roundtrip(self):
        """Model should survive save/load cycle."""
        with tempfile.TemporaryDirectory() as tmpdir:
            model_dir = Path(tmpdir) / "model"
            model_dir.mkdir()

            # Create dummy beta files
            beta1 = model_dir / "beta_site_A_l0.mif"
            beta2 = model_dir / "beta_site_B_l0.mif"
            beta1.touch()
            beta2.touch()

            original = RISHGLMResult(
                site_names=["A", "B"],
                covariate_names=["age"],
                orders=[0, 2],
                beta_paths={
                    "0_site_A": str(beta1),
                    "0_site_B": str(beta2),
                },
                scale_map_paths={},
                reference_site="A",
                n_subjects=10,
                n_per_site={"A": 5, "B": 5},
                cov_means={"age": 35.0},
                cov_stds={"age": 10.0},
                design_columns=["site_A", "site_B", "age"],
                mask_path=None,
                output_dir=str(model_dir),
            )

            json_path = str(model_dir / "model.json")
            save_rish_glm_model(original, json_path)

            loaded = load_rish_glm_model(json_path)

            assert loaded.site_names == original.site_names
            assert loaded.covariate_names == original.covariate_names
            assert loaded.orders == original.orders
            assert loaded.reference_site == original.reference_site
            assert loaded.n_subjects == original.n_subjects
            assert loaded.n_per_site == original.n_per_site
            assert loaded.cov_means == original.cov_means
            assert loaded.cov_stds == original.cov_stds
            assert loaded.design_columns == original.design_columns

    def test_json_format(self):
        """Saved JSON should be valid and readable."""
        with tempfile.TemporaryDirectory() as tmpdir:
            model = RISHGLMResult(
                site_names=["ref", "tar"],
                orders=[0],
                reference_site="ref",
                n_subjects=5,
            )
            json_path = Path(tmpdir) / "model.json"
            save_rish_glm_model(model, str(json_path))

            with open(json_path) as f:
                data = json.load(f)

            assert data["reference_site"] == "ref"
            assert data["site_names"] == ["ref", "tar"]


class TestEdgeCases:
    """Test edge cases in design matrix and fitting."""

    def test_single_subject_per_site(self):
        """Single subject per site should still produce a design matrix."""
        labels = ["A", "B"]
        design, names, _, _, _ = build_rish_glm_design(labels)

        assert design.shape == (2, 2)
        np.testing.assert_array_equal(design, np.eye(2))

    def test_unbalanced_sites(self):
        """Unbalanced site sizes should work correctly."""
        labels = ["A"] * 20 + ["B"] * 3
        design, _, site_map, _, _ = build_rish_glm_design(labels)

        assert design.shape == (23, 2)
        assert design[:, site_map["A"]].sum() == 20
        assert design[:, site_map["B"]].sum() == 3

    def test_many_covariates(self):
        """Design matrix with many covariates."""
        labels = ["A", "A", "B", "B"]
        covariates = {f"cov{i}": [float(j + i) for j in range(4)]
                      for i in range(5)}
        design, names, _, _, _ = build_rish_glm_design(labels, covariates)

        assert design.shape == (4, 7)  # 2 sites + 5 covariates
        assert len(names) == 7

    def test_empty_covariates(self):
        """Empty covariates dict should be treated as None."""
        labels = ["A", "B"]
        design1, _, _, _, _ = build_rish_glm_design(labels, covariates={})
        design2, _, _, _, _ = build_rish_glm_design(labels, covariates=None)

        np.testing.assert_array_equal(design1, design2)

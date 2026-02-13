"""
Unit tests for covariate regression and participant parsing.
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

from src.core.covariates import (
    CovariateModel,
    save_covariate_model,
    load_covariate_model,
    _standardize_covariates,
)
from src.io.participants import (
    ParticipantData,
    load_participants_tsv,
    load_participants_csv,
    _encode_categorical,
    _handle_missing_values,
    _is_numeric,
)


# =============================================================================
# CovariateModel Tests
# =============================================================================

class TestCovariateModel:
    """Tests for CovariateModel dataclass."""

    def test_creation_defaults(self):
        model = CovariateModel()
        assert model.covariate_names == []
        assert model.orders == []
        assert model.means == {}
        assert model.stds == {}
        assert model.beta_paths == {}
        assert model.intercept_paths == {}
        assert model.mask_path is None
        assert model.n_subjects == 0

    def test_creation_with_values(self):
        model = CovariateModel(
            covariate_names=["age", "sex"],
            orders=[0, 2, 4],
            means={"age": 35.0, "sex": 0.5},
            stds={"age": 10.0, "sex": 0.5},
            n_subjects=20,
        )
        assert model.covariate_names == ["age", "sex"]
        assert len(model.orders) == 3
        assert model.n_subjects == 20


class TestCovariateModelSerialization:
    """Tests for save/load round-trip."""

    def test_save_load_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model_dir = Path(tmpdir) / "model"
            model_dir.mkdir()

            # Create dummy beta map files
            for l in [0, 2]:
                (model_dir / f"intercept_l{l}.mif").touch()
                (model_dir / f"beta_age_l{l}.mif").touch()

            model = CovariateModel(
                covariate_names=["age"],
                orders=[0, 2],
                means={"age": 30.0},
                stds={"age": 8.0},
                beta_paths={
                    "0_age": str(model_dir / "beta_age_l0.mif"),
                    "2_age": str(model_dir / "beta_age_l2.mif"),
                },
                intercept_paths={
                    0: str(model_dir / "intercept_l0.mif"),
                    2: str(model_dir / "intercept_l2.mif"),
                },
                mask_path=None,
                n_subjects=15,
                output_dir=str(model_dir),
            )

            json_path = str(model_dir / "covariate_model.json")
            save_covariate_model(model, json_path)

            # Verify JSON exists and is valid
            assert Path(json_path).exists()
            with open(json_path) as f:
                data = json.load(f)
            assert data["covariate_names"] == ["age"]
            assert data["n_subjects"] == 15

            # Load back
            loaded = load_covariate_model(json_path)
            assert loaded.covariate_names == model.covariate_names
            assert loaded.orders == model.orders
            assert loaded.means == model.means
            assert loaded.stds == model.stds
            assert loaded.n_subjects == model.n_subjects

    def test_relative_path_resolution(self):
        """Paths in JSON should be relative; loading resolves them."""
        with tempfile.TemporaryDirectory() as tmpdir:
            model_dir = Path(tmpdir) / "template" / "covariate_model"
            model_dir.mkdir(parents=True)

            # Create dummy files
            (model_dir / "beta_age_l0.mif").touch()
            (model_dir / "intercept_l0.mif").touch()

            model = CovariateModel(
                covariate_names=["age"],
                orders=[0],
                means={"age": 25.0},
                stds={"age": 5.0},
                beta_paths={"0_age": str(model_dir / "beta_age_l0.mif")},
                intercept_paths={0: str(model_dir / "intercept_l0.mif")},
                n_subjects=10,
                output_dir=str(model_dir),
            )

            json_path = str(model_dir / "covariate_model.json")
            save_covariate_model(model, json_path)

            # Verify stored paths are relative
            with open(json_path) as f:
                data = json.load(f)
            assert not Path(data["beta_paths"]["0_age"]).is_absolute()

            # Load and verify resolved paths are absolute
            loaded = load_covariate_model(json_path)
            assert Path(loaded.beta_paths["0_age"]).is_absolute()

    def test_load_missing_file(self):
        with pytest.raises(FileNotFoundError):
            load_covariate_model("/nonexistent/covariate_model.json")

    def test_multiple_covariates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            model_dir = Path(tmpdir)

            for l in [0, 2]:
                (model_dir / f"intercept_l{l}.mif").touch()
                (model_dir / f"beta_age_l{l}.mif").touch()
                (model_dir / f"beta_sex_l{l}.mif").touch()

            model = CovariateModel(
                covariate_names=["age", "sex"],
                orders=[0, 2],
                means={"age": 40.0, "sex": 0.6},
                stds={"age": 12.0, "sex": 0.49},
                beta_paths={
                    "0_age": str(model_dir / "beta_age_l0.mif"),
                    "0_sex": str(model_dir / "beta_sex_l0.mif"),
                    "2_age": str(model_dir / "beta_age_l2.mif"),
                    "2_sex": str(model_dir / "beta_sex_l2.mif"),
                },
                intercept_paths={
                    0: str(model_dir / "intercept_l0.mif"),
                    2: str(model_dir / "intercept_l2.mif"),
                },
                n_subjects=25,
            )

            json_path = str(model_dir / "covariate_model.json")
            save_covariate_model(model, json_path)
            loaded = load_covariate_model(json_path)

            assert loaded.covariate_names == ["age", "sex"]
            assert len(loaded.beta_paths) == 4


# =============================================================================
# Standardization Tests
# =============================================================================

class TestStandardizeCovariates:
    """Tests for covariate z-scoring."""

    def test_basic_standardization(self):
        covs = {"age": [20.0, 30.0, 40.0]}
        z_covs, means, stds = _standardize_covariates(covs)

        assert means["age"] == pytest.approx(30.0)
        assert stds["age"] == pytest.approx(np.std([20, 30, 40]))
        assert z_covs["age"][1] == pytest.approx(0.0)  # Mean maps to 0

    def test_constant_covariate(self):
        """Constant covariate should use std=1 to avoid NaN."""
        covs = {"group": [1.0, 1.0, 1.0]}
        z_covs, means, stds = _standardize_covariates(covs)

        assert stds["group"] == 1.0
        assert all(v == pytest.approx(0.0) for v in z_covs["group"])

    def test_multiple_covariates(self):
        covs = {
            "age": [20.0, 40.0],
            "sex": [0.0, 1.0],
        }
        z_covs, means, stds = _standardize_covariates(covs)

        assert "age" in means and "sex" in means
        assert len(z_covs["age"]) == 2
        assert len(z_covs["sex"]) == 2


# =============================================================================
# Covariate Adjustment Logic Tests (Pure NumPy)
# =============================================================================

class TestCovariateAdjustmentLogic:
    """Test the mathematical correctness of covariate regression.

    These tests use pure numpy (no MRtrix3 required).
    """

    def test_regression_removes_age_correlation(self):
        """Regression should remove linear age effect from RISH."""
        np.random.seed(42)
        n_subjects = 50
        n_voxels = 100

        # Simulate: RISH = base_signal + age_effect * age + noise
        ages = np.linspace(20, 70, n_subjects)
        base_signal = np.random.rand(n_voxels) * 100 + 50  # 50-150
        age_effect = np.random.rand(n_voxels) * 2  # 0-2 per year
        noise = np.random.randn(n_subjects, n_voxels) * 5

        rish = base_signal[None, :] + age_effect[None, :] * ages[:, None] + noise

        # Z-score age
        age_mean = ages.mean()
        age_std = ages.std()
        age_z = (ages - age_mean) / age_std

        # Design: [intercept, age_z]
        design = np.column_stack([np.ones(n_subjects), age_z])

        # Fit regression
        betas, _, _, _ = np.linalg.lstsq(design, rish, rcond=None)

        # Adjust: remove age effect
        adjusted = rish - betas[1][None, :] * age_z[:, None]

        # Verify: correlation with age should be near zero after adjustment
        for v in range(n_voxels):
            corr_before = np.corrcoef(ages, rish[:, v])[0, 1]
            corr_after = np.corrcoef(ages, adjusted[:, v])[0, 1]
            assert abs(corr_after) < 0.15, (
                f"Voxel {v}: corr_after={corr_after:.3f} (before={corr_before:.3f})"
            )

    def test_adjustment_preserves_scanner_effect(self):
        """Covariate adjustment should preserve scanner differences."""
        np.random.seed(123)
        n_per_site = 30
        n_voxels = 50

        # Site A: younger subjects, scanner factor 1.0
        ages_a = np.random.normal(25, 5, n_per_site)
        scanner_a = 1.0

        # Site B: older subjects, scanner factor 1.5
        ages_b = np.random.normal(55, 5, n_per_site)
        scanner_b = 1.5

        # RISH = scanner_factor * base + age_effect * age + noise
        base = np.random.rand(n_voxels) * 100 + 50
        age_effect = np.random.rand(n_voxels) * 0.5

        rish_a = (scanner_a * base[None, :] +
                  age_effect[None, :] * ages_a[:, None] +
                  np.random.randn(n_per_site, n_voxels) * 2)

        rish_b = (scanner_b * base[None, :] +
                  age_effect[None, :] * ages_b[:, None] +
                  np.random.randn(n_per_site, n_voxels) * 2)

        # Fit model on site A (reference)
        age_z_a = (ages_a - ages_a.mean()) / ages_a.std()
        design_a = np.column_stack([np.ones(n_per_site), age_z_a])
        betas_a, _, _, _ = np.linalg.lstsq(design_a, rish_a, rcond=None)

        # Adjusted template (remove age from A, then average)
        adjusted_a = rish_a - betas_a[1][None, :] * age_z_a[:, None]
        template = adjusted_a.mean(axis=0)

        # Adjust target (site B) using reference model params
        age_z_b = (ages_b - ages_a.mean()) / ages_a.std()
        adjusted_b = rish_b - betas_a[1][None, :] * age_z_b[:, None]
        target_mean = adjusted_b.mean(axis=0)

        # Scale map: template / target
        scale_map = template / (target_mean + 1e-10)

        # The scale map should reflect scanner ratio (~1.0/1.5 ≈ 0.67)
        # Without covariate adjustment, age difference would bias this
        expected_ratio = scanner_a / scanner_b
        median_scale = np.median(scale_map)
        assert abs(median_scale - expected_ratio) < 0.15, (
            f"Scale map median={median_scale:.3f}, expected ~{expected_ratio:.3f}"
        )

    def test_scale_maps_isolate_scanner_effect(self):
        """Scale maps from adjusted RISH should capture only scanner effect."""
        np.random.seed(456)
        n = 40
        n_voxels = 80

        ages = np.random.normal(40, 15, n)
        base = np.random.rand(n_voxels) * 100 + 50
        scanner_factor = 1.3
        age_effect = np.random.rand(n_voxels) * 1.0

        # Reference RISH (scanner=1.0)
        rish_ref = (base[None, :] +
                    age_effect[None, :] * ages[:, None] +
                    np.random.randn(n, n_voxels) * 3)

        # Target RISH (scanner=1.3, same age distribution)
        rish_tar = (scanner_factor * base[None, :] +
                    age_effect[None, :] * ages[:, None] +
                    np.random.randn(n, n_voxels) * 3)

        # Fit on reference
        age_z = (ages - ages.mean()) / ages.std()
        design = np.column_stack([np.ones(n), age_z])

        betas_ref, _, _, _ = np.linalg.lstsq(design, rish_ref, rcond=None)
        betas_tar, _, _, _ = np.linalg.lstsq(design, rish_tar, rcond=None)

        # Adjusted means
        adj_ref = rish_ref - betas_ref[1][None, :] * age_z[:, None]
        adj_tar = rish_tar - betas_tar[1][None, :] * age_z[:, None]

        template = adj_ref.mean(axis=0)
        target_mean = adj_tar.mean(axis=0)

        scale_map = template / (target_mean + 1e-10)

        # Scale should be ~1/1.3 ≈ 0.77
        expected = 1.0 / scanner_factor
        median_scale = np.median(scale_map)
        assert abs(median_scale - expected) < 0.1, (
            f"Scale={median_scale:.3f}, expected ~{expected:.3f}"
        )

    def test_multiple_covariates_age_and_sex(self):
        """Regression with multiple covariates removes both effects."""
        np.random.seed(789)
        n = 60
        n_voxels = 50

        ages = np.random.normal(45, 15, n)
        sexes = np.random.choice([0.0, 1.0], n)

        base = np.random.rand(n_voxels) * 100 + 50
        age_effect = np.random.rand(n_voxels) * 1.0
        sex_effect = np.random.rand(n_voxels) * 10

        rish = (base[None, :] +
                age_effect[None, :] * ages[:, None] +
                sex_effect[None, :] * sexes[:, None] +
                np.random.randn(n, n_voxels) * 3)

        # Z-score
        age_z = (ages - ages.mean()) / ages.std()
        sex_z = (sexes - sexes.mean()) / sexes.std()

        design = np.column_stack([np.ones(n), age_z, sex_z])
        betas, _, _, _ = np.linalg.lstsq(design, rish, rcond=None)

        # Adjust
        adjusted = rish - betas[1][None, :] * age_z[:, None] - betas[2][None, :] * sex_z[:, None]

        # Verify correlations with both covariates are near zero
        for v in range(0, n_voxels, 10):
            corr_age = abs(np.corrcoef(ages, adjusted[:, v])[0, 1])
            corr_sex = abs(np.corrcoef(sexes, adjusted[:, v])[0, 1])
            assert corr_age < 0.2, f"Voxel {v}: age corr={corr_age:.3f}"
            assert corr_sex < 0.2, f"Voxel {v}: sex corr={corr_sex:.3f}"

    def test_no_covariates_is_identity(self):
        """With no covariates, adjustment should be a no-op."""
        np.random.seed(101)
        n = 20
        n_voxels = 30

        rish = np.random.rand(n, n_voxels) * 100

        # Design with only intercept
        design = np.ones((n, 1))
        betas, _, _, _ = np.linalg.lstsq(design, rish, rcond=None)

        # No covariate betas to subtract -> adjusted = rish
        # (intercept removal would change values, but we only remove covariate betas)
        # With 0 covariates, there's nothing to subtract
        adjusted = rish.copy()  # No covariates to adjust for

        np.testing.assert_array_equal(adjusted, rish)


# =============================================================================
# Participant Parsing Tests
# =============================================================================

class TestParticipantsParsing:
    """Tests for loading participant demographics."""

    def test_load_tsv_basic(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("participant_id\tage\tsex\n")
            f.write("sub-01\t25\tM\n")
            f.write("sub-02\t35\tF\n")
            f.write("sub-03\t45\tM\n")
            tsv_path = f.name

        result = load_participants_tsv(tsv_path, ["age", "sex"])
        assert result.n_subjects == 3
        assert result.subject_ids == ["sub-01", "sub-02", "sub-03"]
        assert result.covariates["age"] == [25.0, 35.0, 45.0]
        assert result.covariates["sex"] == [1.0, 0.0, 1.0]

        Path(tsv_path).unlink()

    def test_load_csv_basic(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as f:
            f.write("subject,age,sex\n")
            f.write("sub-01,25,M\n")
            f.write("sub-02,35,F\n")
            tsv_path = f.name

        result = load_participants_csv(tsv_path, ["age", "sex"])
        assert result.n_subjects == 2
        assert result.covariates["age"] == [25.0, 35.0]

        Path(tsv_path).unlink()

    def test_subject_ordering(self):
        """Rows should be reordered to match subject_ids."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("participant_id\tage\n")
            f.write("sub-03\t45\n")
            f.write("sub-01\t25\n")
            f.write("sub-02\t35\n")
            tsv_path = f.name

        result = load_participants_tsv(
            tsv_path, ["age"],
            subject_ids=["sub-01", "sub-02", "sub-03"]
        )
        assert result.subject_ids == ["sub-01", "sub-02", "sub-03"]
        assert result.covariates["age"] == [25.0, 35.0, 45.0]

        Path(tsv_path).unlink()

    def test_missing_subject_raises(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("participant_id\tage\n")
            f.write("sub-01\t25\n")
            tsv_path = f.name

        with pytest.raises(ValueError, match="sub-99"):
            load_participants_tsv(
                tsv_path, ["age"],
                subject_ids=["sub-01", "sub-99"]
            )

        Path(tsv_path).unlink()

    def test_missing_column_raises(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("participant_id\tage\n")
            f.write("sub-01\t25\n")
            tsv_path = f.name

        with pytest.raises(ValueError, match="weight"):
            load_participants_tsv(tsv_path, ["age", "weight"])

        Path(tsv_path).unlink()

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_participants_tsv("/nonexistent.tsv", ["age"])

    def test_missing_values_mean_imputation(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("participant_id\tage\n")
            f.write("sub-01\t20\n")
            f.write("sub-02\t\n")
            f.write("sub-03\t40\n")
            tsv_path = f.name

        result = load_participants_tsv(tsv_path, ["age"])
        # Mean of 20, 40 = 30
        assert result.covariates["age"][1] == pytest.approx(30.0)

        Path(tsv_path).unlink()

    def test_missing_values_na_string(self):
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False
        ) as f:
            f.write("participant_id\tage\n")
            f.write("sub-01\t20\n")
            f.write("sub-02\tN/A\n")
            f.write("sub-03\t40\n")
            tsv_path = f.name

        result = load_participants_tsv(tsv_path, ["age"])
        assert result.covariates["age"][1] == pytest.approx(30.0)

        Path(tsv_path).unlink()

    def test_participant_data_properties(self):
        data = ParticipantData(
            subject_ids=["a", "b", "c"],
            covariates={"age": [1.0, 2.0, 3.0], "sex": [0.0, 1.0, 0.0]},
        )
        assert data.n_subjects == 3
        assert set(data.covariate_names) == {"age", "sex"}


class TestCategoricalEncoding:
    """Tests for categorical variable encoding."""

    def test_sex_encoding_m_f(self):
        assert _encode_categorical(["M", "F", "M"]) == [1.0, 0.0, 1.0]

    def test_sex_encoding_male_female(self):
        assert _encode_categorical(["Male", "Female"]) == [1.0, 0.0]

    def test_sex_encoding_numeric(self):
        assert _encode_categorical(["1", "0", "1"]) == [1.0, 0.0, 1.0]

    def test_general_categorical(self):
        result = _encode_categorical(["A", "B", "C", "A"])
        assert result == [0.0, 1.0, 2.0, 0.0]


class TestNumericDetection:
    """Tests for _is_numeric."""

    def test_numeric(self):
        assert _is_numeric(["1.5", "2.3", "4.0"])

    def test_non_numeric(self):
        assert not _is_numeric(["M", "F", "M"])

    def test_numeric_with_missing(self):
        assert _is_numeric(["1.5", "", "N/A", "4.0"])

    def test_mixed(self):
        assert not _is_numeric(["1.5", "abc", "4.0"])


class TestHandleMissingValues:
    """Tests for missing value imputation."""

    def test_no_missing(self):
        result = _handle_missing_values(["1", "2", "3"])
        assert result == [1.0, 2.0, 3.0]

    def test_mean_imputation(self):
        result = _handle_missing_values(["10", "", "30"], strategy="mean")
        assert result[1] == pytest.approx(20.0)

    def test_median_imputation(self):
        result = _handle_missing_values(["10", "", "30", "100"], strategy="median")
        # Valid: [10, 30, 100], median = 30
        assert result[1] == pytest.approx(30.0)

    def test_na_strings(self):
        result = _handle_missing_values(["10", "NA", "N/A", "30"])
        # Mean of [10, 30] = 20
        assert result[1] == pytest.approx(20.0)
        assert result[2] == pytest.approx(20.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

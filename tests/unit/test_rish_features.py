"""
Unit tests for RISH feature extraction.
"""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest

from src.core.rish_features import get_sh_indices, SHInfo, extract_rish_features


def mrtrix_available():
    """Check if MRtrix3 is available."""
    try:
        result = subprocess.run(
            ["mrinfo", "--version"],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def nibabel_available():
    """Check if nibabel is available."""
    try:
        import nibabel
        return True
    except ImportError:
        return False


class TestSHIndices:
    """Tests for SH index calculation."""
    
    def test_lmax_0(self):
        """Test lmax=0 (just DC component)."""
        info = get_sh_indices(0)
        assert info.lmax == 0
        assert info.n_volumes == 1
        assert info.volume_indices == {0: (0, 1)}
        assert info.n_coeffs_per_order == {0: 1}
    
    def test_lmax_2(self):
        """Test lmax=2."""
        info = get_sh_indices(2)
        assert info.lmax == 2
        assert info.n_volumes == 6  # 1 + 5
        assert info.volume_indices == {0: (0, 1), 2: (1, 6)}
        assert info.n_coeffs_per_order == {0: 1, 2: 5}
    
    def test_lmax_4(self):
        """Test lmax=4."""
        info = get_sh_indices(4)
        assert info.lmax == 4
        assert info.n_volumes == 15  # 1 + 5 + 9
        assert info.volume_indices == {
            0: (0, 1),
            2: (1, 6),
            4: (6, 15)
        }
    
    def test_lmax_6(self):
        """Test lmax=6."""
        info = get_sh_indices(6)
        assert info.lmax == 6
        assert info.n_volumes == 28  # 1 + 5 + 9 + 13
    
    def test_lmax_8(self):
        """Test lmax=8 (common case)."""
        info = get_sh_indices(8)
        assert info.lmax == 8
        assert info.n_volumes == 45  # 1 + 5 + 9 + 13 + 17
        assert info.volume_indices == {
            0: (0, 1),
            2: (1, 6),
            4: (6, 15),
            6: (15, 28),
            8: (28, 45)
        }
    
    def test_odd_lmax_raises(self):
        """Test that odd lmax raises error."""
        with pytest.raises(ValueError, match="must be even"):
            get_sh_indices(3)
        
        with pytest.raises(ValueError, match="must be even"):
            get_sh_indices(7)
    
    def test_coeffs_per_order(self):
        """Test coefficient counts: 2l+1 for each order."""
        info = get_sh_indices(8)
        
        assert info.n_coeffs_per_order[0] == 1   # 2*0+1
        assert info.n_coeffs_per_order[2] == 5   # 2*2+1
        assert info.n_coeffs_per_order[4] == 9   # 2*4+1
        assert info.n_coeffs_per_order[6] == 13  # 2*6+1
        assert info.n_coeffs_per_order[8] == 17  # 2*8+1


# Path to existing test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_output"


def _test_data_available():
    """Check if test data is available (helper function)."""
    fod_path = TEST_DATA_DIR / "sub-002" / "ses-1" / "wm_fod.mif"
    return fod_path.exists()


@pytest.mark.skipif(
    not mrtrix_available() or not _test_data_available(),
    reason="MRtrix3 or test data not available"
)
class TestRISHExtraction:
    """Tests for RISH feature extraction (requires MRtrix3 and test data)."""

    @pytest.fixture
    def real_sh_data(self):
        """Use existing FOD data from test_output."""
        sh_path = TEST_DATA_DIR / "sub-002" / "ses-1" / "wm_fod.mif"
        lmax = 8  # 45 volumes = lmax 8
        return sh_path, lmax

    @pytest.fixture
    def real_mask(self):
        """Use existing mask from test_output."""
        return TEST_DATA_DIR / "sub-002" / "ses-1" / "mask.mif"

    @pytest.fixture
    def existing_rish(self):
        """Reference to already extracted RISH features."""
        rish_dir = TEST_DATA_DIR / "sub-002" / "ses-1" / "rish"
        return {
            0: rish_dir / "rish_l0.mif",
            2: rish_dir / "rish_l2.mif",
            4: rish_dir / "rish_l4.mif",
            6: rish_dir / "rish_l6.mif",
            8: rish_dir / "rish_l8.mif",
        }

    def test_extract_basic(self, real_sh_data, tmp_path):
        """Test basic RISH extraction without mask."""
        sh_path, lmax = real_sh_data
        output_dir = tmp_path / "rish_output"

        rish = extract_rish_features(
            str(sh_path),
            str(output_dir),
            lmax=lmax,
            mask=None,
            n_threads=1
        )

        # Check outputs exist for each order (lmax=8 -> l=0,2,4,6,8)
        assert 0 in rish
        assert 2 in rish
        assert 4 in rish
        assert 6 in rish
        assert 8 in rish

        # Check files exist
        for l, path in rish.items():
            assert Path(path).exists(), f"RISH l={l} file not found"

        # Verify output dimensions (should be 3D or 4D with size 1 on 4th dim)
        result = subprocess.run(
            ["mrinfo", "-size", rish[0]],
            capture_output=True, text=True, check=True
        )
        dims = [int(x) for x in result.stdout.strip().split()]
        # MRtrix may keep 4th dimension with size 1
        assert len(dims) >= 3, f"Expected at least 3D output, got {len(dims)}D"
        if len(dims) == 4:
            assert dims[3] == 1, f"Expected 4th dim to be 1, got {dims[3]}"

    def test_extract_with_mask(self, real_sh_data, real_mask, tmp_path):
        """Test RISH extraction with brain mask."""
        sh_path, lmax = real_sh_data
        output_dir = tmp_path / "rish_output_masked"

        rish = extract_rish_features(
            str(sh_path),
            str(output_dir),
            lmax=lmax,
            mask=str(real_mask),
            n_threads=1
        )

        # Check outputs exist
        assert len(rish) == 5  # l=0, 2, 4, 6, 8

        # Check files exist
        for l, path in rish.items():
            assert Path(path).exists()

    def test_rish_values_positive(self, existing_rish):
        """Test that RISH values are non-negative (they are norms)."""
        # Check existing RISH files
        for l, path in existing_rish.items():
            result = subprocess.run(
                ["mrstats", str(path), "-output", "min"],
                capture_output=True, text=True, check=True
            )
            min_val = float(result.stdout.strip())
            assert min_val >= 0.0, f"RISH l={l} has negative values (min={min_val})"

    def test_extract_lmax_autodetect(self, real_sh_data, tmp_path):
        """Test RISH extraction with auto-detected lmax."""
        sh_path, expected_lmax = real_sh_data
        output_dir = tmp_path / "rish_auto"

        # Don't specify lmax - should auto-detect from 45 volumes -> lmax=8
        rish = extract_rish_features(
            str(sh_path),
            str(output_dir),
            lmax=None,  # Auto-detect
            n_threads=1
        )

        # Should detect lmax=8 from 45 volumes
        assert 8 in rish
        assert 10 not in rish  # Should not have l=10

    def test_rish_consistency(self, real_sh_data, real_mask, existing_rish, tmp_path):
        """Test that newly extracted RISH matches existing RISH."""
        sh_path, lmax = real_sh_data
        output_dir = tmp_path / "rish_consistency"

        # Extract fresh RISH
        rish = extract_rish_features(
            str(sh_path),
            str(output_dir),
            lmax=lmax,
            mask=str(real_mask),
            n_threads=1
        )

        # Compare with existing (should be very similar)
        for l in [0, 2, 4, 6, 8]:
            # Get stats from both
            result_new = subprocess.run(
                ["mrstats", rish[l], "-output", "mean"],
                capture_output=True, text=True, check=True
            )
            mean_new = float(result_new.stdout.strip())

            result_old = subprocess.run(
                ["mrstats", str(existing_rish[l]), "-output", "mean"],
                capture_output=True, text=True, check=True
            )
            mean_old = float(result_old.stdout.strip())

            # Should be identical (or very close)
            assert abs(mean_new - mean_old) < 0.01 * max(mean_new, mean_old, 1e-6), \
                f"RISH l={l} mismatch: new={mean_new}, old={mean_old}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

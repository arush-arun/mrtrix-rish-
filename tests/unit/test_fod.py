"""
Unit tests for FOD computation module.
"""

import pytest
from src.core.fod import ShellType, ShellInfo


class TestShellInfo:
    """Tests for ShellInfo dataclass."""
    
    def test_single_shell(self):
        """Test single-shell detection."""
        info = ShellInfo(
            shell_type=ShellType.SINGLE_SHELL,
            bvalues=[1000],
            shell_counts={1000: 64},
            n_b0=6,
            n_total=70
        )
        
        assert info.n_shells == 1
        assert not info.is_multi_shell
        assert info.shell_type == ShellType.SINGLE_SHELL
    
    def test_multi_shell(self):
        """Test multi-shell detection."""
        info = ShellInfo(
            shell_type=ShellType.MULTI_SHELL,
            bvalues=[1000, 2000, 3000],
            shell_counts={1000: 64, 2000: 64, 3000: 64},
            n_b0=18,
            n_total=210
        )
        
        assert info.n_shells == 3
        assert info.is_multi_shell
        assert info.shell_type == ShellType.MULTI_SHELL
    
    def test_hcp_style(self):
        """Test HCP-style acquisition (3 shells)."""
        info = ShellInfo(
            shell_type=ShellType.MULTI_SHELL,
            bvalues=[1000, 2000, 3000],
            shell_counts={1000: 90, 2000: 90, 3000: 90},
            n_b0=18,
            n_total=288
        )
        
        assert info.n_shells == 3
        assert info.is_multi_shell


class TestShellDetection:
    """Tests for shell detection (requires MRtrix3)."""
    
    @pytest.mark.skip(reason="Requires MRtrix3 and test data")
    def test_detect_single_shell(self):
        """Test detecting single-shell data."""
        pass
    
    @pytest.mark.skip(reason="Requires MRtrix3 and test data")
    def test_detect_multi_shell(self):
        """Test detecting multi-shell data."""
        pass


class TestFODComputation:
    """Tests for FOD computation (requires MRtrix3)."""
    
    @pytest.mark.skip(reason="Requires MRtrix3 and test data")
    def test_csd_single_shell(self):
        """Test CSD on single-shell data."""
        pass
    
    @pytest.mark.skip(reason="Requires MRtrix3 and test data")
    def test_msmt_csd_multi_shell(self):
        """Test MSMT-CSD on multi-shell data."""
        pass


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

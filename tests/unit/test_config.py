"""
Unit tests for configuration handling.
"""

import pytest
import tempfile
from pathlib import Path

from src.io.config_io import (
    HarmonizationConfig,
    RegistrationConfig,
    QCConfig,
    PipelineConfig,
    load_config,
    save_config,
    get_default_config
)


class TestDefaultConfig:
    """Tests for default configuration."""
    
    def test_default_harmonization(self):
        """Test default harmonization config."""
        config = HarmonizationConfig()
        assert config.lmax == 8
        assert config.smoothing_fwhm == 3.0
        assert config.clip_range == (0.5, 2.0)
        assert config.n_threads == 4
    
    def test_default_registration(self):
        """Test default registration config."""
        config = RegistrationConfig()
        assert config.method == "mrtrix"
        assert config.template == "population"
    
    def test_default_qc(self):
        """Test default QC config."""
        config = QCConfig()
        assert config.generate_report == True
        assert "fa_diff" in config.metrics
    
    def test_get_default_config(self):
        """Test getting full default config."""
        config = get_default_config()
        assert isinstance(config, PipelineConfig)
        assert isinstance(config.harmonization, HarmonizationConfig)


class TestConfigIO:
    """Tests for config load/save."""
    
    def test_save_and_load(self):
        """Test round-trip save and load."""
        config = get_default_config()
        config.harmonization.lmax = 6
        config.harmonization.smoothing_fwhm = 5.0
        
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "test_config.yaml"
            save_config(config, str(config_path))
            
            loaded = load_config(str(config_path))
            
            assert loaded.harmonization.lmax == 6
            assert loaded.harmonization.smoothing_fwhm == 5.0
            assert loaded.harmonization.clip_range == (0.5, 2.0)
    
    def test_load_missing_file(self):
        """Test loading non-existent file."""
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/config.yaml")
    
    def test_clip_range_conversion(self):
        """Test that clip_range converts from list to tuple."""
        config = get_default_config()
        
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "test_config.yaml"
            save_config(config, str(config_path))
            
            loaded = load_config(str(config_path))
            
            # Should be tuple, not list
            assert isinstance(loaded.harmonization.clip_range, tuple)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

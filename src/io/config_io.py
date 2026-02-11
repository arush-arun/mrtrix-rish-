"""
Configuration I/O Module

Load and save YAML configuration files.
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass, field, asdict


@dataclass
class HarmonizationConfig:
    """Configuration for RISH harmonization."""
    
    # SH parameters
    lmax: int = 8
    
    # Scale map parameters
    mask_threshold: float = 0.1
    smoothing_fwhm: float = 3.0
    clip_range: tuple = (0.5, 2.0)
    
    # Processing
    n_threads: int = 4
    
    # Output
    save_intermediate: bool = False
    compress_output: bool = True


@dataclass 
class RegistrationConfig:
    """Configuration for registration."""
    
    method: str = "mrtrix"  # or "ants"
    template: str = "population"  # or path to template
    linear_type: str = "affine"
    nonlinear: bool = True


@dataclass
class QCConfig:
    """Configuration for QC reports."""
    
    generate_report: bool = True
    metrics: list = field(default_factory=lambda: [
        "fa_diff", "md_diff", "angular_correlation"
    ])
    save_figures: bool = True


@dataclass
class PipelineConfig:
    """Full pipeline configuration."""
    
    harmonization: HarmonizationConfig = field(default_factory=HarmonizationConfig)
    registration: RegistrationConfig = field(default_factory=RegistrationConfig)
    qc: QCConfig = field(default_factory=QCConfig)


def load_config(config_path: str) -> PipelineConfig:
    """
    Load configuration from YAML file.
    
    Parameters
    ----------
    config_path : str
        Path to YAML config file
        
    Returns
    -------
    PipelineConfig
        Loaded configuration
    """
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path) as f:
        data = yaml.safe_load(f)
    
    # Parse into dataclasses
    harm_data = data.get("harmonization", {})
    reg_data = data.get("registration", {})
    qc_data = data.get("qc", {})
    
    # Handle tuple conversion for clip_range
    if "clip_range" in harm_data and isinstance(harm_data["clip_range"], list):
        harm_data["clip_range"] = tuple(harm_data["clip_range"])
    
    config = PipelineConfig(
        harmonization=HarmonizationConfig(**harm_data),
        registration=RegistrationConfig(**reg_data),
        qc=QCConfig(**qc_data)
    )
    
    return config


def save_config(config: PipelineConfig, output_path: str) -> None:
    """
    Save configuration to YAML file.
    
    Parameters
    ----------
    config : PipelineConfig
        Configuration to save
    output_path : str
        Output path
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert dataclasses to dict
    data = {
        "harmonization": asdict(config.harmonization),
        "registration": asdict(config.registration),
        "qc": asdict(config.qc)
    }
    
    # Convert tuple to list for YAML
    data["harmonization"]["clip_range"] = list(data["harmonization"]["clip_range"])
    
    with open(output_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def get_default_config() -> PipelineConfig:
    """Get default configuration."""
    return PipelineConfig()

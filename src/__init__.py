"""
mrtrix-rish: MRtrix3-native RISH Harmonization Pipeline

A reproducible tool for multi-site diffusion MRI harmonization using
Rotational Invariant Spherical Harmonics (RISH).
"""

__version__ = "0.1.0"
__author__ = "Arush"

from .core.rish_features import extract_rish_features
from .core.scale_maps import compute_scale_maps
from .core.harmonize import harmonize_sh, RISHHarmonizer

__all__ = [
    "extract_rish_features",
    "compute_scale_maps", 
    "harmonize_sh",
    "RISHHarmonizer",
]

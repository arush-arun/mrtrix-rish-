"""
QC and reporting modules.

Includes:
- metrics: Basic QC metrics (FA/MD differences)
- reports: HTML report generation
- glm: Voxel-wise General Linear Model (following MRtrix3 design)
- site_effects: Site effect statistical testing with permutation inference
"""

from .metrics import compute_qc_metrics
from .reports import generate_html_report
from .glm import (
    Hypothesis,
    Partition,
    TestFixedHomoscedastic,
    TestFixedHeteroscedastic,
    GLMResult,
    TestOutput,
    create_design_matrix,
    create_site_contrast,
    check_design,
)
from .site_effects import (
    test_site_effect,
    compare_site_effects,
    SiteEffectResult,
    Shuffler,
    fdr_correction,
    compute_partial_eta_squared,
)

__all__ = [
    # Metrics
    "compute_qc_metrics",
    "generate_html_report",
    # GLM
    "Hypothesis",
    "Partition",
    "TestFixedHomoscedastic",
    "TestFixedHeteroscedastic",
    "GLMResult",
    "TestOutput",
    "create_design_matrix",
    "create_site_contrast",
    "check_design",
    # Site effects
    "test_site_effect",
    "compare_site_effects",
    "SiteEffectResult",
    "Shuffler",
    "fdr_correction",
    "compute_partial_eta_squared",
]

# mrtrix-rish Design Document

## 1. Overview

### 1.1 Problem Statement
Multi-site diffusion MRI data exhibits scanner-dependent variability that confounds group analyses. Current harmonization methods either:
- Operate at the scalar level (FA, MD) — losing angular information
- Require FSL/ANTs dependencies — complicating MRtrix3 workflows

### 1.2 Solution
Implement RISH harmonization natively in MRtrix3, operating on spherical harmonic coefficients to preserve fiber orientation information while removing site effects.

### 1.3 Design Principles
1. **Modularity** — Each component is independent and testable
2. **Reproducibility** — Containerized, version-controlled, BIDS-compatible
3. **Transparency** — Clear logging, QC outputs, provenance tracking
4. **Performance** — Leverage MRtrix3's multi-threading

---

## 2. Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        mrtrix-rish                              │
├─────────────────────────────────────────────────────────────────┤
│  CLI Layer (cli/)                                               │
│  ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐              │
│  │harmonize│ │ create  │ │  apply  │ │   qc    │              │
│  │         │ │template │ │         │ │         │              │
│  └────┬────┘ └────┬────┘ └────┬────┘ └────┬────┘              │
├───────┼──────────┼──────────┼──────────┼────────────────────────┤
│  Core Layer (core/)                                             │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐               │
│  │ rish_features│ │ scale_maps  │ │ harmonize   │               │
│  │  .py        │ │  .py        │ │  .py        │               │
│  └─────────────┘ └─────────────┘ └─────────────┘               │
├─────────────────────────────────────────────────────────────────┤
│  I/O Layer (io/)                                                │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐               │
│  │ mrtrix_io   │ │ bids_io     │ │ config_io   │               │
│  └─────────────┘ └─────────────┘ └─────────────┘               │
├─────────────────────────────────────────────────────────────────┤
│  QC Layer (qc/)                                                 │
│  ┌─────────────┐ ┌─────────────┐ ┌─────────────┐               │
│  │ reports     │ │ visualize   │ │ metrics     │               │
│  └─────────────┘ └─────────────┘ └─────────────┘               │
├─────────────────────────────────────────────────────────────────┤
│  External: MRtrix3 binaries (amp2sh, sh2amp, mrcalc, mrmath)   │
└─────────────────────────────────────────────────────────────────┘
```

---

## 3. Core Components

### 3.0 FOD Computation (`core/fod.py`)

**Purpose:** Automatically compute FODs using the appropriate algorithm based on shell structure.

**Shell Detection:**
```python
shell_info = detect_shells("dwi.mif")
# Returns: ShellInfo(shell_type, bvalues, shell_counts, n_b0, n_total)
```

**Algorithm Selection:**
| Shell Structure | Algorithm | Response Functions |
|-----------------|-----------|-------------------|
| Single-shell (1 b-value) | CSD | WM only |
| Multi-shell (2+ b-values) | MSMT-CSD | WM, GM, CSF |

**MRtrix3 Commands:**
- Shell detection: `mrinfo -shell_bvalues`, `mrinfo -shell_sizes`
- Response estimation: `dwi2response dhollander` (multi-shell) or `dwi2response tournier` (single-shell)
- FOD computation: `dwi2fod csd` or `dwi2fod msmt_csd`

**Example:**
```python
from mrtrix_rish.core import compute_fod

# Auto-detects shells and picks CSD or MSMT-CSD
result = compute_fod(
    dwi="subject_dwi.mif",
    mask="brain_mask.mif",
    output_dir="fod_output/"
)
# result["wm_fod"] → WM FOD image
# result["algorithm"] → "csd" or "msmt_csd"
```

### 3.1 RISH Feature Extraction (`core/rish_features.py`)

**Purpose:** Extract rotationally invariant features (θₗ) from SH coefficients.

**Math:**
```
θₗ = √(Σₘ |cₗₘ|²)   for m ∈ [-l, +l]
```

**MRtrix3 Implementation:**
- Input: SH coefficient image (from `amp2sh`)
- For each order l (0, 2, 4, 6, 8...):
  - Extract relevant volumes using `mrconvert -coord 3 <indices>`
  - Square coefficients: `mrcalc <input> -sq`
  - Sum across m: `mrmath <input> sum -axis 3`
  - Square root: `mrcalc <input> -sqrt`
- Output: θ₀, θ₂, θ₄, θ₆... images

**SH Volume Indexing:**
| Order (l) | Volumes | Count |
|-----------|---------|-------|
| 0         | 0       | 1     |
| 2         | 1-5     | 5     |
| 4         | 6-14    | 9     |
| 6         | 15-27   | 13    |
| 8         | 28-44   | 17    |

### 3.2 Template Creation (`core/template.py`)

**Purpose:** Create reference RISH features in common space.

**Steps:**
1. Register all subjects to common space (population template or MNI)
2. Compute RISH features for each subject
3. Average RISH features across reference site subjects
4. Output: Template θₗ images

### 3.3 Scale Map Computation (`core/scale_maps.py`)

**Purpose:** Compute per-voxel, per-order scaling factors.

**Math:**
```
scale_l(x) = θₗ_ref(x) / θₗ_target(x)
```

**Considerations:**
- Handle division by zero (mask low-signal regions)
- Smooth scale maps to avoid noise amplification
- Clip extreme values

### 3.4 Harmonization (`core/harmonize.py`)

**Purpose:** Apply scaling to SH coefficients.

**Steps:**
1. Load target SH coefficients
2. For each order l:
   - Multiply all m coefficients by scale_l
3. Reconstruct harmonized signal (optional, via `sh2amp`)

---

## 4. Data Flow

```
INPUT                    PROCESSING                     OUTPUT
─────                    ──────────                     ──────

Reference Site           
┌──────────────┐        ┌───────────────┐
│ DWI + mask   │───────▶│   amp2sh      │
│ (N subjects) │        └───────┬───────┘
└──────────────┘                │
                                ▼
                        ┌───────────────┐
                        │ Extract RISH  │──────▶ θₗ_ref (template)
                        │   features    │
                        └───────────────┘

Target Site                     
┌──────────────┐        ┌───────────────┐
│ DWI + mask   │───────▶│   amp2sh      │
│ (M subjects) │        └───────┬───────┘
└──────────────┘                │
                                ▼
                        ┌───────────────┐
                        │ Extract RISH  │──────▶ θₗ_target
                        │   features    │
                        └───────────────┘
                                │
                                ▼
                        ┌───────────────┐
                        │ Compute scale │──────▶ scale_l = θₗ_ref / θₗ_target
                        │     maps      │
                        └───────────────┘
                                │
                                ▼
                        ┌───────────────┐        ┌──────────────┐
                        │    Apply      │───────▶│ Harmonized   │
                        │   scaling     │        │ SH / DWI     │
                        └───────────────┘        └──────────────┘
                                │
                                ▼
                        ┌───────────────┐        ┌──────────────┐
                        │   Generate    │───────▶│  QC Report   │
                        │      QC       │        │   (HTML)     │
                        └───────────────┘        └──────────────┘
```

---

## 5. File Structure

```
mrtrix-rish/
├── src/
│   ├── __init__.py
│   ├── core/
│   │   ├── __init__.py
│   │   ├── rish_features.py    # θₗ extraction
│   │   ├── template.py         # Template creation
│   │   ├── scale_maps.py       # Scale computation
│   │   └── harmonize.py        # Apply harmonization
│   ├── io/
│   │   ├── __init__.py
│   │   ├── mrtrix_io.py        # MRtrix3 file handling
│   │   ├── bids_io.py          # BIDS dataset handling
│   │   └── config_io.py        # YAML/JSON config
│   ├── qc/
│   │   ├── __init__.py
│   │   ├── metrics.py          # Quantitative metrics
│   │   ├── visualize.py        # Plots and figures
│   │   └── reports.py          # HTML report generation
│   └── cli/
│       ├── __init__.py
│       ├── main.py             # Entry point
│       ├── harmonize.py        # harmonize subcommand
│       ├── template.py         # create-template subcommand
│       └── qc.py               # qc subcommand
├── tests/
│   ├── unit/                   # Unit tests per module
│   ├── integration/            # Full pipeline tests
│   └── data/                   # Test fixtures
├── docs/
│   ├── design/                 # This document
│   └── api/                    # Auto-generated API docs
├── scripts/
│   ├── install_mrtrix.sh       # MRtrix3 installation helper
│   └── validate_env.sh         # Check dependencies
├── config/
│   ├── default.yaml            # Default parameters
│   └── example_harmonize.yaml  # Example config
├── examples/
│   ├── minimal/                # Minimal working example
│   └── hcp_multisite/          # HCP-style multi-site example
├── docker/
│   ├── Dockerfile
│   └── Singularity.def
├── pyproject.toml
├── LICENSE
└── README.md
```

---

## 6. Configuration

### 6.1 Example Config File

```yaml
# config/harmonize.yaml
harmonization:
  lmax: 8                    # Maximum SH order
  mask_threshold: 0.1        # Minimum signal for masking
  scale_smoothing_fwhm: 3.0  # Gaussian smoothing (mm)
  scale_clip_range: [0.5, 2.0]  # Min/max scaling factors

registration:
  method: "mrtrix"           # or "ants"
  template: "population"     # or path to external template

qc:
  generate_report: true
  metrics: ["fa_diff", "md_diff", "angular_correlation"]
  
output:
  save_intermediate: true
  compress: true
```

---

## 7. QC Metrics

| Metric | Description | Target |
|--------|-------------|--------|
| FA difference | Mean |ΔFA| before/after | < 0.02 |
| MD difference | Mean |ΔMD| before/after | < 5% |
| Angular correlation | FOD similarity (ACC) | > 0.9 |
| Scale map range | % voxels with extreme scaling | < 5% |
| Site effect (GLM) | Residual site variance | p > 0.05 |

---

## 8. API Design

### 8.1 Python API

```python
from mrtrix_rish import RISHHarmonizer

# Initialize
harmonizer = RISHHarmonizer(
    lmax=8,
    n_threads=4
)

# Create template from reference site
harmonizer.create_template(
    reference_dwi_list=["site_a/sub-01/dwi.mif", ...],
    reference_mask_list=["site_a/sub-01/mask.mif", ...],
    output_template="template/"
)

# Harmonize target site
harmonizer.harmonize(
    target_dwi="site_b/sub-01/dwi.mif",
    target_mask="site_b/sub-01/mask.mif",
    template="template/",
    output="harmonized/sub-01/"
)

# Generate QC report
harmonizer.generate_qc_report(
    original="site_b/sub-01/dwi.mif",
    harmonized="harmonized/sub-01/dwi.mif",
    output="qc/sub-01_report.html"
)
```

### 8.2 CLI

```bash
# Create template
mrtrix-rish create-template \
    --reference-list reference_subjects.txt \
    --output template/ \
    --lmax 8

# Harmonize
mrtrix-rish harmonize \
    --target site_b/sub-01/dwi.mif \
    --mask site_b/sub-01/mask.mif \
    --template template/ \
    --output harmonized/

# QC only
mrtrix-rish qc \
    --original site_b/sub-01/dwi.mif \
    --harmonized harmonized/sub-01/dwi.mif \
    --output qc/
```

---

## 9. Testing Strategy

### 9.1 Unit Tests
- RISH feature extraction on synthetic SH data
- Scale map computation with known inputs
- I/O functions with mock data

### 9.2 Integration Tests
- Full pipeline on HCP test-retest data
- Comparison with pnlbwh implementation
- BIDS dataset processing

### 9.3 Validation
- Replicate Cetin Karayumak 2019 results
- Statistical test for site effect removal

---

## 10. Milestones

| Phase | Deliverable | Duration |
|-------|-------------|----------|
| 1 | Core RISH extraction + unit tests | 1 week |
| 2 | Template creation + scale maps | 1 week |
| 3 | Harmonization + CLI | 1 week |
| 4 | QC reports + documentation | 1 week |
| 5 | Docker/Singularity + validation | 1 week |
| 6 | Paper draft + release | 2 weeks |

---

## 11. References

1. Cetin Karayumak et al. (2019). Retrospective harmonization of multi-site diffusion MRI data. NeuroImage.
2. Mirzaalian et al. (2018). Multi-site harmonization of diffusion MRI data in a registration framework. Brain Imaging Behav.
3. De Luca et al. (2025). Cross-site harmonization of diffusion MRI data without matched training subjects. MRM.
4. Tournier et al. (2019). MRtrix3: A fast, flexible and open software framework. NeuroImage.

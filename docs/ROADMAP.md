# mrtrix-rish Development Roadmap

> Version: 1.0
> Last Updated: 2026-02-02
> Focus Areas: Covariate Support, Regional QC, Statistical Site Effect Testing, Tractography Validation, Tutorials

---

## Executive Summary

This roadmap outlines the development plan to address critical gaps in mrtrix-rish, transforming it from an alpha-stage tool to a production-ready solution for the neuroimaging community. The focus areas are:

| Priority | Feature | Impact | Effort |
|----------|---------|--------|--------|
| 1 | Covariate Support | High | Medium |
| 2 | Statistical Site Effect Test | High | Low |
| 3 | Regional QC | Medium | Medium |
| 4 | Tractography Validation | High | High |
| 5 | Tutorials & Documentation | High | Medium |

---

## Phase 1: Statistical Site Effect Testing (Weeks 1-2)

### Rationale
Without statistical validation, users cannot objectively assess whether harmonization successfully removed site effects. This is the foundation for all other QC features.

### 1.1 GLM-Based Site Effect Test

**File:** `src/qc/site_effects.py` (new)

**Implementation:**

```python
# Core function signature
def test_site_effect(
    metric_images: Dict[str, List[str]],  # {site_id: [fa_paths...]}
    mask: str,
    covariates: Optional[pd.DataFrame] = None,
    test_type: str = "anova",  # "anova", "permutation", "combat_check"
    n_permutations: int = 5000,
    output_dir: Optional[str] = None
) -> SiteEffectResult
```

**Statistical Methods:**

| Method | Use Case | Output |
|--------|----------|--------|
| ANOVA F-test | Quick assessment | F-statistic, p-value per voxel |
| Permutation test | Robust inference | Cluster-corrected p-values |
| ICC (Intraclass Correlation) | Reliability assessment | ICC(1,1) maps |
| Effect size (Cohen's d) | Magnitude assessment | d-maps between sites |

**Algorithm:**

```
1. Load FA/MD images for all subjects grouped by site
2. Fit voxel-wise GLM: metric ~ site + covariates
3. Extract F-statistic for site factor
4. Apply multiple comparison correction (FDR or permutation)
5. Generate summary:
   - % voxels with significant site effect (p < 0.05)
   - Mean effect size across brain
   - Cluster-level statistics
```

**Success Criteria:**
- Before harmonization: Significant site effect expected
- After harmonization: <5% voxels with p < 0.05 (at FDR-corrected threshold)
- Effect size reduction: >80% reduction in mean Cohen's d

### 1.2 Integration with CLI

**New CLI command:**

```bash
mrtrix-rish site-effect \
    --site-list sites.csv \        # CSV with columns: subject, site, fa_path, [age, sex]
    --mask template_mask.mif \
    --test permutation \
    --n-permutations 5000 \
    --covariates age,sex \
    --output site_effect_report/
```

### 1.3 Deliverables

| Deliverable | Description |
|-------------|-------------|
| `src/qc/site_effects.py` | Core statistical testing module |
| `src/qc/glm.py` | Voxel-wise GLM implementation |
| `tests/unit/test_site_effects.py` | Unit tests |
| `docs/tutorials/site_effect_testing.md` | User guide |

---

## Phase 2: Covariate Support (Weeks 3-5)

### Rationale
Biological variation (age, sex) is confounded with site effects in cross-sectional studies. Without covariate modeling, harmonization may remove real biological signal.

### 2.1 Architecture Decision

**Option A: Pre-harmonization Covariate Regression (Recommended)**
- Remove age/sex effects before RISH computation
- Re-add after harmonization
- Preserves biological signal

**Option B: Covariate-Adjusted RISH Templates**
- Create age/sex-stratified templates
- More complex, requires larger samples

**Decision:** Implement Option A with Option B as future extension.

### 2.2 Covariate Module

**File:** `src/core/covariates.py` (new)

```python
@dataclass
class CovariateModel:
    """Stores fitted covariate model for preservation."""
    coefficients: Dict[str, np.ndarray]  # {covariate: beta_map}
    residual_image: str
    mean_image: str

def fit_covariates(
    images: List[str],
    covariates: pd.DataFrame,
    mask: str,
    output_dir: str,
    variables: List[str] = ["age", "sex"]
) -> CovariateModel:
    """
    Fit voxel-wise regression: image ~ age + sex + intercept
    Store coefficients for later restoration.
    """

def remove_covariates(
    image: str,
    model: CovariateModel,
    subject_covariates: Dict[str, float],
    output: str
) -> str:
    """Remove covariate effects: residual = image - (age*beta_age + sex*beta_sex)"""

def restore_covariates(
    harmonized_image: str,
    model: CovariateModel,
    subject_covariates: Dict[str, float],
    output: str
) -> str:
    """Restore covariate effects to harmonized image."""
```

### 2.3 Updated Workflow

```
Original Workflow:
DWI → SH → RISH → Scale → Harmonize

New Workflow with Covariates:
DWI → SH → Remove Covariates → RISH → Scale → Harmonize → Restore Covariates
         ↓                                                    ↑
    Store betas ──────────────────────────────────────────────┘
```

### 2.4 Configuration Extension

**config/default.yaml additions:**

```yaml
covariates:
  enabled: false
  variables:
    - name: age
      type: continuous
      standardize: true
    - name: sex
      type: categorical
      reference: "M"
  model_file: null  # Path to pre-fitted model
  fit_from_reference: true  # Fit model from reference site
```

### 2.5 BIDS Integration

**Participant.tsv parsing:**

```python
def load_bids_covariates(bids_dir: str) -> pd.DataFrame:
    """
    Load covariates from participants.tsv

    Expected columns:
    - participant_id
    - age (years)
    - sex (M/F)
    - site (optional)
    """
```

### 2.6 Deliverables

| Deliverable | Description |
|-------------|-------------|
| `src/core/covariates.py` | Covariate regression module |
| `src/io/participants.py` | BIDS participants.tsv parser |
| `tests/unit/test_covariates.py` | Unit tests |
| CLI flag `--covariates` | Enable covariate correction |
| `docs/tutorials/covariate_correction.md` | User guide |

---

## Phase 3: Regional QC (Weeks 6-8)

### Rationale
Whole-brain metrics mask regional failures. Atlas-based QC reveals if specific structures (e.g., corpus callosum) are poorly harmonized.

### 3.1 Atlas Integration

**Supported Atlases:**

| Atlas | Regions | Use Case |
|-------|---------|----------|
| JHU ICBM-DTI-81 | 48 WM tracts | Standard DTI analysis |
| JHU White Matter | 20 WM regions | Tract-based statistics |
| Desikan-Killiany | 68 cortical | GM reference |
| AAL3 | 170 regions | Comprehensive coverage |
| Custom | User-defined | Study-specific ROIs |

**File:** `src/qc/regional.py` (new)

```python
@dataclass
class RegionalQCResult:
    atlas_name: str
    regions: Dict[str, RegionMetrics]
    summary_table: pd.DataFrame
    flagged_regions: List[str]  # Regions failing QC

@dataclass
class RegionMetrics:
    name: str
    n_voxels: int
    fa_pre: float
    fa_post: float
    fa_diff: float
    md_pre: float
    md_post: float
    md_diff: float
    site_effect_p: float
    pass_qc: bool

def compute_regional_qc(
    pre_harmonization: str,
    post_harmonization: str,
    atlas: str,
    mask: str,
    output_dir: str
) -> RegionalQCResult:
    """Compute per-region harmonization metrics."""
```

### 3.2 Regional Statistics Pipeline

```
1. Register atlas to subject space (or subject to template)
2. For each region:
   a. Extract voxels within region mask
   b. Compute pre/post FA, MD statistics
   c. Compute effect size of change
   d. Test for residual site effect within region
3. Flag regions with:
   - |ΔFA| > 0.03
   - |ΔMD| > 10%
   - Residual site effect p < 0.05
4. Generate regional report
```

### 3.3 Visualization

**Regional QC Report Elements:**

1. **Glass brain plot** - Color-coded regions by QC status
2. **Table** - Per-region statistics with pass/fail indicators
3. **Radar chart** - Multi-metric comparison across key tracts
4. **Histogram** - Distribution of changes per region

### 3.4 CLI Integration

```bash
mrtrix-rish qc-regional \
    --pre original_fa.mif \
    --post harmonized_fa.mif \
    --atlas JHU-ICBM-DTI-81 \
    --mask brain_mask.mif \
    --output regional_qc/
```

### 3.5 Deliverables

| Deliverable | Description |
|-------------|-------------|
| `src/qc/regional.py` | Regional QC computation |
| `src/qc/atlases.py` | Atlas management and registration |
| `data/atlases/` | Bundled atlas files (JHU, AAL3) |
| `src/qc/regional_report.py` | HTML report with visualizations |
| `tests/unit/test_regional_qc.py` | Unit tests |

---

## Phase 4: Tractography Validation (Weeks 9-13)

### Rationale
RISH harmonization claims to preserve fiber orientation. This must be validated by demonstrating that tractography results are consistent across sites after harmonization.

### 4.1 Validation Framework

**File:** `src/validation/tractography.py` (new)

```python
@dataclass
class TractographyValidationResult:
    tract_name: str
    pre_harmonization: TractMetrics
    post_harmonization: TractMetrics
    site_variability_reduction: float
    angular_error: float

@dataclass
class TractMetrics:
    streamline_count: int
    mean_length: float
    volume: float
    fa_along_tract: np.ndarray
    endpoint_density: np.ndarray
```

### 4.2 Validation Pipeline

```
Phase A: Reference Tract Generation
1. Run tractography on reference site subjects
2. Create tract probability maps (tract density imaging)
3. Define "gold standard" tract templates

Phase B: Cross-Site Validation
For each site (before and after harmonization):
1. Run identical tractography protocol
2. Compute tract metrics:
   - Dice overlap with reference tracts
   - Streamline count variability (CV across sites)
   - Mean FA along tract profile
   - Angular correlation of FOD peaks
3. Statistical comparison:
   - Site effect on tract metrics (ANOVA)
   - ICC for tract reproducibility
```

### 4.3 Tract-Specific Metrics

**Tracts to Validate:**

| Tract | Rationale |
|-------|-----------|
| Corpus Callosum (genu, body, splenium) | High SNR, easy to track |
| Corticospinal Tract | Clinical relevance |
| Arcuate Fasciculus | Language studies |
| Cingulum | Limbic system |
| Uncinate Fasciculus | Frontal-temporal |
| Superior Longitudinal Fasciculus | Association fibers |

### 4.4 Implementation

```python
def validate_tractography(
    reference_fods: List[str],
    target_fods_pre: List[str],
    target_fods_post: List[str],
    tract_definitions: Dict[str, TractDefinition],
    output_dir: str,
    tractography_params: TractographyParams = None
) -> TractographyValidationReport:
    """
    Full tractography validation pipeline.

    Steps:
    1. Generate tracks using tckgen
    2. Filter with tckedit using tract ROIs
    3. Compute tract density images
    4. Extract along-tract profiles
    5. Compare pre/post harmonization
    """
```

### 4.5 MRtrix3 Commands Used

```bash
# FOD-based tractography
tckgen fod.mif tracks.tck -seed_image seed.mif -include roi.mif -select 10000

# Tract filtering
tckedit tracks.tck filtered.tck -include roi1.mif -include roi2.mif

# Tract density imaging
tckmap tracks.tck tdi.mif -template fa.mif

# Along-tract sampling
tcksample tracks.tck fa.mif fa_along_tract.csv

# Track statistics
tckstats tracks.tck
```

### 4.6 Deliverables

| Deliverable | Description |
|-------------|-------------|
| `src/validation/tractography.py` | Tractography validation module |
| `src/validation/tract_profiles.py` | Along-tract analysis |
| `src/validation/tract_definitions.py` | Standard tract ROI definitions |
| `data/tract_rois/` | Tract seed/target ROIs |
| `tests/integration/test_tractography.py` | Integration tests |
| `docs/tutorials/tractography_validation.md` | User guide |

---

## Phase 5: Tutorials & Documentation (Weeks 14-16)

### Rationale
Adoption depends on usability. Comprehensive tutorials lower the barrier to entry and demonstrate best practices.

### 5.1 Tutorial Structure

```
docs/tutorials/
├── 01_getting_started/
│   ├── installation.md
│   ├── quick_start.md
│   └── example_data.md
├── 02_basic_harmonization/
│   ├── single_site_pair.md
│   ├── multi_site.md
│   └── bids_workflow.md
├── 03_advanced_features/
│   ├── covariate_correction.md
│   ├── custom_templates.md
│   └── shell_mismatch.md
├── 04_quality_control/
│   ├── site_effect_testing.md
│   ├── regional_qc.md
│   └── tractography_validation.md
├── 05_troubleshooting/
│   ├── common_errors.md
│   ├── faq.md
│   └── performance_tuning.md
└── 06_best_practices/
    ├── study_design.md
    ├── reference_site_selection.md
    └── reporting_guidelines.md
```

### 5.2 Tutorial 1: Getting Started

**docs/tutorials/01_getting_started/quick_start.md**

```markdown
# Quick Start Guide

## Prerequisites
- MRtrix3 >= 3.0.4
- Python >= 3.9
- Preprocessed DWI data (denoised, motion/distortion corrected)

## Installation

### Option 1: pip
pip install mrtrix-rish

### Option 2: Docker
docker pull arush/mrtrix-rish:latest

### Option 3: From source
git clone https://github.com/arush/mrtrix-rish
cd mrtrix-rish
pip install -e ".[qc]"

## Verify Installation
mrtrix-rish --version
mrtrix-rish --help

## Download Example Data
mrtrix-rish download-examples --output examples/

## Run Your First Harmonization
# 1. Create template from reference site
mrtrix-rish create-template \
    --reference-list examples/reference_subjects.txt \
    --output examples/template/

# 2. Harmonize target site
mrtrix-rish harmonize \
    --target examples/target/sub-01_dwi.mif \
    --mask examples/target/sub-01_mask.mif \
    --template examples/template/ \
    --output examples/harmonized/

# 3. Generate QC report
mrtrix-rish qc \
    --original examples/target/sub-01_dwi.mif \
    --harmonized examples/harmonized/sh_harmonized.mif \
    --output examples/qc_report.html
```

### 5.3 Tutorial 2: BIDS Workflow

**docs/tutorials/02_basic_harmonization/bids_workflow.md**

```markdown
# Processing a BIDS Dataset

## Dataset Structure
your_study/
├── participants.tsv          # Must include 'site' column
├── sub-01/
│   └── ses-01/
│       └── dwi/
│           ├── sub-01_ses-01_dwi.nii.gz
│           ├── sub-01_ses-01_dwi.bval
│           └── sub-01_ses-01_dwi.bvec
├── sub-02/
...

## Step 1: Identify Sites
mrtrix-rish bids-list /path/to/bids --group-by site

## Step 2: Select Reference Site
# Choose site with:
# - Largest sample size
# - Best data quality
# - Most representative scanner

## Step 3: Run BIDS Harmonization
mrtrix-rish bids /path/to/bids \
    --output /path/to/derivatives/mrtrix-rish \
    --reference-site "SiteA" \
    --covariates age,sex \
    --lmax 8 \
    --n-threads 8

## Step 4: Validate Results
mrtrix-rish site-effect \
    --derivatives /path/to/derivatives/mrtrix-rish \
    --test permutation \
    --output site_effect_report/
```

### 5.4 Tutorial 3: Covariate Correction

**docs/tutorials/03_advanced_features/covariate_correction.md**

```markdown
# Covariate Correction for Age and Sex

## Why Use Covariates?
Site effects are confounded with demographics. If Site A has older subjects:
- Without covariates: Age effect removed as "site effect"
- With covariates: Age effect preserved, only scanner effect removed

## Prepare participants.tsv
participant_id  age     sex     site
sub-01          25      M       SiteA
sub-02          67      F       SiteA
sub-03          45      M       SiteB
...

## Run with Covariates
mrtrix-rish bids /path/to/bids \
    --output derivatives/ \
    --covariates age,sex \
    --reference-site SiteA

## Verify Covariate Preservation
# After harmonization, age effects should still be present:
mrtrix-rish analyze-covariates \
    --derivatives derivatives/ \
    --variable age \
    --output covariate_analysis/

# Expected: Significant age-FA correlation preserved
```

### 5.5 Tutorial 4: Site Effect Testing

**docs/tutorials/04_quality_control/site_effect_testing.md**

```markdown
# Statistical Testing for Site Effects

## Concept
Before harmonization: Significant site effect (expected)
After harmonization: No significant site effect (goal)

## Run Site Effect Test

### Before Harmonization
mrtrix-rish site-effect \
    --input-list pre_harmonization.csv \
    --mask template_mask.mif \
    --test permutation \
    --n-permutations 5000 \
    --output site_effect_pre/

### After Harmonization
mrtrix-rish site-effect \
    --input-list post_harmonization.csv \
    --mask template_mask.mif \
    --test permutation \
    --n-permutations 5000 \
    --output site_effect_post/

## Interpret Results

### site_effect_pre/summary.json
{
    "test": "permutation_anova",
    "n_sites": 3,
    "n_subjects": 150,
    "significant_voxels_percent": 45.2,  # HIGH = site effect present
    "mean_effect_size": 0.82,
    "conclusion": "SIGNIFICANT_SITE_EFFECT"
}

### site_effect_post/summary.json
{
    "test": "permutation_anova",
    "n_sites": 3,
    "n_subjects": 150,
    "significant_voxels_percent": 3.1,   # LOW = harmonization successful
    "mean_effect_size": 0.15,
    "effect_size_reduction": 0.82,       # 82% reduction
    "conclusion": "NO_SIGNIFICANT_SITE_EFFECT"
}

## QC Thresholds
| Metric | Pass | Warning | Fail |
|--------|------|---------|------|
| % significant voxels | <5% | 5-15% | >15% |
| Effect size reduction | >70% | 50-70% | <50% |
```

### 5.6 Tutorial 5: Regional QC

**docs/tutorials/04_quality_control/regional_qc.md**

```markdown
# Regional Quality Control with Atlases

## Available Atlases
- JHU-ICBM-DTI-81: 48 white matter tract labels
- JHU-WhiteMatter: 20 white matter regions
- AAL3: 170 brain regions

## Run Regional QC
mrtrix-rish qc-regional \
    --pre fa_original.mif \
    --post fa_harmonized.mif \
    --atlas JHU-ICBM-DTI-81 \
    --output regional_qc/

## Output Files
regional_qc/
├── regional_metrics.csv      # Per-region statistics
├── flagged_regions.txt       # Regions failing QC
├── regional_report.html      # Visual report
└── figures/
    ├── glass_brain.png
    ├── radar_chart.png
    └── region_histograms.png

## Interpret Results

### regional_metrics.csv
region,n_voxels,fa_pre,fa_post,fa_diff,site_p,qc_pass
Genu_CC,1523,0.72,0.71,0.01,0.42,PASS
Body_CC,2841,0.68,0.67,0.01,0.38,PASS
CST_L,892,0.61,0.58,0.03,0.02,FAIL  # <-- Flagged
...

## Troubleshooting Flagged Regions
1. Check registration quality in that region
2. Verify mask coverage
3. Consider region-specific harmonization
```

### 5.7 Tutorial 6: Tractography Validation

**docs/tutorials/04_quality_control/tractography_validation.md**

```markdown
# Validating Tractography After Harmonization

## Purpose
RISH harmonization preserves fiber orientations. Validate by comparing tractography results across sites.

## Run Validation
mrtrix-rish validate-tractography \
    --reference-fods reference_site/*.mif \
    --target-fods-pre target_site_pre/*.mif \
    --target-fods-post target_site_post/*.mif \
    --tracts CST,CC,AF \
    --output tractography_validation/

## Key Metrics

### Dice Overlap with Reference
| Tract | Pre-Harm | Post-Harm | Improvement |
|-------|----------|-----------|-------------|
| CST_L | 0.72 | 0.89 | +24% |
| CC | 0.81 | 0.91 | +12% |
| AF_L | 0.68 | 0.85 | +25% |

### Cross-Site Variability (CV)
| Tract | Pre-Harm CV | Post-Harm CV | Reduction |
|-------|-------------|--------------|-----------|
| CST_L | 18% | 6% | 67% |
| CC | 12% | 4% | 67% |
| AF_L | 22% | 8% | 64% |

## Interpretation
- Dice > 0.85: Excellent tract agreement
- CV < 10%: Acceptable cross-site variability
- CV reduction > 50%: Harmonization effective
```

### 5.8 Example Dataset

**Provide downloadable example data:**

```
examples/
├── README.md
├── download.sh               # Script to fetch data
├── site_A/                   # Reference site (5 subjects)
│   ├── sub-01/
│   │   ├── dwi.mif
│   │   └── mask.mif
│   └── ...
├── site_B/                   # Target site (5 subjects)
│   └── ...
├── participants.tsv          # Demographics
├── expected_outputs/         # For validation
│   ├── template/
│   ├── harmonized/
│   └── qc_reports/
└── run_example.sh           # End-to-end script
```

### 5.9 Deliverables

| Deliverable | Description |
|-------------|-------------|
| `docs/tutorials/` | Complete tutorial set (6 sections) |
| `examples/` | Downloadable example dataset |
| `examples/run_example.sh` | End-to-end example script |
| `docs/API.md` | Python API reference |
| `docs/CLI.md` | CLI reference with all flags |
| Video tutorials | YouTube walkthroughs (optional) |

---

## Implementation Timeline

```
Week    Phase                           Milestone
──────────────────────────────────────────────────────────────
1-2     Statistical Site Effect Test    GLM module complete
                                        CLI integration
                                        Unit tests passing
──────────────────────────────────────────────────────────────
3-5     Covariate Support               Covariate regression module
                                        BIDS integration
                                        Workflow updated
──────────────────────────────────────────────────────────────
6-8     Regional QC                     Atlas integration
                                        Regional metrics
                                        HTML reports
──────────────────────────────────────────────────────────────
9-13    Tractography Validation         Tract generation pipeline
                                        Validation metrics
                                        Integration tests
──────────────────────────────────────────────────────────────
14-16   Tutorials & Documentation       All tutorials written
                                        Example data available
                                        API documentation
──────────────────────────────────────────────────────────────
```

---

## File Structure After Implementation

```
src/
├── core/
│   ├── harmonize.py          # Existing
│   ├── rish_features.py      # Existing
│   ├── scale_maps.py         # Existing
│   ├── covariates.py         # NEW - Phase 2
│   └── ...
├── qc/
│   ├── metrics.py            # Existing (enhanced)
│   ├── reports.py            # Existing (enhanced)
│   ├── site_effects.py       # NEW - Phase 1
│   ├── glm.py                # NEW - Phase 1
│   ├── regional.py           # NEW - Phase 3
│   └── atlases.py            # NEW - Phase 3
├── validation/
│   ├── tractography.py       # NEW - Phase 4
│   ├── tract_profiles.py     # NEW - Phase 4
│   └── tract_definitions.py  # NEW - Phase 4
├── io/
│   ├── participants.py       # NEW - Phase 2
│   └── ...
└── cli/
    └── main.py               # Updated with new commands

data/
├── atlases/
│   ├── JHU-ICBM-DTI-81.nii.gz
│   ├── JHU-WhiteMatter.nii.gz
│   └── labels/
└── tract_rois/
    ├── CST/
    ├── CC/
    └── AF/

docs/
├── tutorials/                # NEW - Phase 5
│   ├── 01_getting_started/
│   ├── 02_basic_harmonization/
│   ├── 03_advanced_features/
│   ├── 04_quality_control/
│   ├── 05_troubleshooting/
│   └── 06_best_practices/
├── API.md                    # NEW - Phase 5
└── CLI.md                    # NEW - Phase 5

examples/                     # NEW - Phase 5
├── site_A/
├── site_B/
├── participants.tsv
└── run_example.sh
```

---

## Dependencies to Add

```toml
# pyproject.toml additions

[project.optional-dependencies]
stats = [
    "scipy>=1.9.0",           # Statistical tests
    "statsmodels>=0.14.0",    # GLM, mixed effects
    "pingouin>=0.5.0",        # ICC computation
]
regional = [
    "nibabel>=4.0",           # Atlas I/O
    "nilearn>=0.10.0",        # Atlas utilities
]
validation = [
    "pandas>=2.0",            # Data handling
    "seaborn>=0.12",          # Visualization
]
full = [
    "mrtrix-rish[qc,stats,regional,validation]"
]
```

---

## Testing Strategy

### Unit Tests

| Module | Test File | Coverage Target |
|--------|-----------|-----------------|
| site_effects.py | test_site_effects.py | 90% |
| covariates.py | test_covariates.py | 90% |
| regional.py | test_regional.py | 85% |
| tractography.py | test_tractography.py | 80% |

### Integration Tests

| Test | Description |
|------|-------------|
| test_covariate_workflow.py | End-to-end with covariates |
| test_full_qc_pipeline.py | All QC modules together |
| test_bids_with_covariates.py | BIDS + covariates |

### Validation Tests

| Test | Data Required |
|------|---------------|
| test_site_effect_removal.py | Multi-site synthetic data |
| test_covariate_preservation.py | Age-varying synthetic data |
| test_tractography_consistency.py | FOD test data |

---

## Success Criteria

### Phase 1: Site Effect Testing
- [ ] GLM-based test correctly identifies known site effects
- [ ] Permutation test matches parametric results
- [ ] <5% false positive rate on harmonized data

### Phase 2: Covariate Support
- [ ] Age-FA correlation preserved after harmonization
- [ ] Sex differences maintained
- [ ] No interaction between covariates and harmonization

### Phase 3: Regional QC
- [ ] All JHU-81 regions correctly segmented
- [ ] Flagged regions match manual inspection
- [ ] HTML report renders correctly

### Phase 4: Tractography Validation
- [ ] CST Dice > 0.85 after harmonization
- [ ] Cross-site CV < 10%
- [ ] Angular correlation > 0.9

### Phase 5: Tutorials
- [ ] New user can complete quick start in <30 minutes
- [ ] All code examples execute without errors
- [ ] Example dataset downloads correctly

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| Covariate model overfitting | Cross-validation, regularization |
| Atlas registration failures | Multiple registration attempts, manual QC flag |
| Tractography variability | Standardized parameters, multiple seeds |
| Large atlas files | Lazy loading, optional download |
| Tutorial code rot | CI testing of all examples |

---

## Future Extensions (Out of Scope)

These are noted for future roadmaps:

1. **Deep learning harmonization** - Neural network-based approach
2. **Longitudinal mode** - Within-subject harmonization
3. **Traveling phantom validation** - Gold-standard testing
4. **GPU acceleration** - CuPy/JAX backend
5. **Multi-reference templates** - Population average reference
6. **Interactive QC dashboard** - Web-based visualization

---

## Contact & Contributions

- **Issues:** https://github.com/arush/mrtrix-rish/issues
- **Discussions:** https://github.com/arush/mrtrix-rish/discussions
- **Contributing:** See CONTRIBUTING.md

---

*This roadmap is a living document and will be updated as development progresses.*

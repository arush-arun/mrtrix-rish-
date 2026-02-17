# mrtrix-rish

**MRtrix3-native RISH Harmonization Pipeline**

A reproducible, open-source tool for multi-site diffusion MRI harmonization using Rotationally Invariant Spherical Harmonics (RISH). Implements both the classical RISH method (Mirzaalian et al., 2016; Cetin Karayumak et al., 2019) and the RISH-GLM joint modeling approach (De Luca et al., 2025).

## Why?

Multi-site diffusion MRI studies suffer from scanner-induced variability that confounds downstream analyses. RISH harmonization operates at the spherical harmonics level, preserving angular (fiber orientation) information while removing site effects. Existing implementations (e.g., [pnlbwh/dMRIharmonization](https://github.com/pnlbwh/dMRIharmonization))

## Features

- **Two harmonization methods:**
  - **Classical RISH** — two-step template creation + per-subject harmonization (Mirzaalian et al., 2016; Cetin Karayumak et al., 2019), with optional separate covariate adjustment
  - **RISH-GLM** — joint GLM estimation of site effects and covariates in a single model (De Luca et al., 2025), no matched training subjects required
- **FOD-level harmonization** — preserves fiber orientation information
- **Multi-site support** — harmonize across 2+ sites simultaneously (RISH-GLM)
- **Covariate modeling** — account for age, sex, and other confounds
- **Built-in site-effect testing** — GLM with permutation inference to validate harmonization
- **BIDS-compatible** — works with standard neuroimaging data structures

## Prerequisites

All FOD images must be **registered to a common template space** before harmonization. A typical preprocessing pipeline:

1. Compute FODs from DWI data (`dwi2response` + `dwi2fod`)
2. Register all FODs to a reference subject or population template (`mrregister` + `mrtransform`)
3. Create a group brain mask (`mrmath ... min`)

## Methods

### Classical RISH (two-step)

The classical approach (Mirzaalian et al., 2016; Cetin Karayumak et al., 2019) works in two stages:

1. **Create template** from reference site: extract RISH features, optionally adjust for covariates, average across subjects
2. **Harmonize** each target subject: compute per-voxel, per-SH-order scale maps `θ_l(x) = template_l(x) / target_l(x)`, apply to SH coefficients

**Step 1: Create a reference list** — a text file with one registered FOD path per line:

```
# ref_fods.txt
/data/site_A/sub-01/fod_reg.mif
/data/site_A/sub-02/fod_reg.mif
/data/site_A/sub-03/fod_reg.mif
```

Similarly for masks:

```
# ref_masks.txt
/data/site_A/sub-01/mask_reg.mif
/data/site_A/sub-02/mask_reg.mif
/data/site_A/sub-03/mask_reg.mif
```

**Step 2: Create the template and harmonize:**

```bash
# Build template from reference site
mrtrix-rish create-template \
    --reference-list ref_fods.txt \
    --mask-list ref_masks.txt \
    --output template/ \
    --lmax 8

# Harmonize a target subject
mrtrix-rish harmonize \
    --target /data/site_B/sub-01/fod_reg.mif \
    --template template/ \
    --mask /data/site_B/sub-01/mask_reg.mif \
    --output harmonized/sub-01/
```

**With covariate adjustment** (optional):

```bash
# Template with covariate regression
mrtrix-rish create-template \
    --reference-list ref_fods.txt \
    --mask-list ref_masks.txt \
    --participants participants.tsv \
    --covariates age,sex \
    --output template_cov/

# Harmonize with subject-specific covariates
mrtrix-rish harmonize \
    --target /data/site_B/sub-01/fod_reg.mif \
    --template template_cov/ \
    --mask /data/site_B/sub-01/mask_reg.mif \
    --subject-covariates "age=35,sex=M" \
    --output harmonized/sub-01/
```

### RISH-GLM (joint model)

The RISH-GLM approach (De Luca et al., 2025) fits a single GLM per voxel across all subjects from all sites simultaneously:

```
RISH_l(v) = β_R(v) * S_R + β_T(v) * S_T + β_cov(v) * covariates + ε
```

Scale factors are derived directly from the fitted coefficients: `θ_l(v) = β_R(v) / β_T(v)`.

Advantages over the classical approach:
- No matched training subjects required — different subjects at each site
- Multi-site simultaneous harmonization (>2 sites in one command)
- Joint covariate-site estimation in a single model
- Preserves biological variability with unmatched groups

#### Manifest CSV format

The RISH-GLM command takes a single CSV file describing all subjects. **Required columns:**

| Column | Description |
|--------|-------------|
| `subject` | Unique subject identifier |
| `site` | Site/scanner name |
| `fod_path` | Path to registered FOD image (in template space) |

**Optional covariate columns** can be added (specify which to use with `--covariates`):

| Column | Description |
|--------|-------------|
| `age` | Subject age (numeric) |
| `sex` | Subject sex (M/F, encoded automatically) |
| Any other | Additional numeric or categorical covariates |

**Example manifest — 2 sites, no covariates:**

```csv
subject,site,fod_path
sub-01_siteA,SiteA,/data/siteA/sub-01/fod_reg.mif
sub-02_siteA,SiteA,/data/siteA/sub-02/fod_reg.mif
sub-03_siteA,SiteA,/data/siteA/sub-03/fod_reg.mif
sub-01_siteB,SiteB,/data/siteB/sub-01/fod_reg.mif
sub-02_siteB,SiteB,/data/siteB/sub-02/fod_reg.mif
sub-03_siteB,SiteB,/data/siteB/sub-03/fod_reg.mif
```

**Example manifest — 3 sites with covariates:**

```csv
subject,site,fod_path,age,sex
sub-01_GE,GE_MR750,/data/ge/sub-01/fod_reg.mif,32.5,F
sub-02_GE,GE_MR750,/data/ge/sub-02/fod_reg.mif,45.1,M
sub-03_GE,GE_MR750,/data/ge/sub-03/fod_reg.mif,28.3,F
sub-01_Siemens,Siemens_Prisma,/data/siemens/sub-01/fod_reg.mif,38.7,M
sub-02_Siemens,Siemens_Prisma,/data/siemens/sub-02/fod_reg.mif,55.2,F
sub-03_Siemens,Siemens_Prisma,/data/siemens/sub-03/fod_reg.mif,22.9,M
sub-01_Philips,Philips_Achieva,/data/philips/sub-01/fod_reg.mif,41.0,F
sub-02_Philips,Philips_Achieva,/data/philips/sub-02/fod_reg.mif,33.8,M
```

**Note:** Subjects do NOT need to be the same across sites (unlike the classical method). Each site can have different subjects and different sample sizes.

#### Running RISH-GLM

```bash
# Basic: 2 sites, no covariates
mrtrix-rish rish-glm \
    --manifest manifest.csv \
    --reference-site SiteA \
    --mask group_mask.mif \
    --output output/

# With covariates and automatic harmonization of all target subjects
mrtrix-rish rish-glm \
    --manifest manifest.csv \
    --reference-site GE_MR750 \
    --mask group_mask.mif \
    --covariates age,sex \
    --output output/ \
    --harmonize

# Custom smoothing and clipping
mrtrix-rish rish-glm \
    --manifest manifest.csv \
    --reference-site SiteA \
    --mask group_mask.mif \
    --output output/ \
    --smoothing-fwhm 4.0 \
    --clip-range 0.3,3.0 \
    --lmax 6 \
    --harmonize
```

The `--harmonize` flag applies the computed scale maps to all non-reference-site subjects and writes harmonized FODs to `output/harmonized/<site>/<subject>_harmonized.mif`.

### Site-effect testing

Validate harmonization by testing for residual site effects using a voxel-wise GLM with permutation inference:

```bash
# Create a site list CSV (subject, site, image_path columns)
# image_path points to a scalar metric (e.g., RISH l0 feature)
mrtrix-rish site-effect \
    --site-list sites.csv \
    --mask group_mask.mif \
    --output results/ \
    --n-permutations 5000
```

The site list CSV format:

```csv
subject,site,image_path
sub-01_siteA,SiteA,/data/siteA/sub-01/rish_l0.mif
sub-02_siteA,SiteA,/data/siteA/sub-02/rish_l0.mif
sub-01_siteB,SiteB,/data/siteB/sub-01/rish_l0.mif
sub-02_siteB,SiteB,/data/siteB/sub-02/rish_l0.mif
```

## Quick Start

```bash
# Install
pip install -e .

# Classical harmonization
mrtrix-rish create-template -r ref_fods.txt -m ref_masks.txt -o template/
mrtrix-rish harmonize -t target_fod.mif -T template/ -m mask.mif -o output/

# RISH-GLM harmonization (all sites at once)
mrtrix-rish rish-glm -i manifest.csv -r ReferenceSite -m mask.mif -o output/ --harmonize

# Test for site effects
mrtrix-rish site-effect -s site_list.csv -m mask.mif -o results/
```

## Requirements

- MRtrix3 >= 3.0.4
- Python >= 3.9
- NumPy, NiBabel, SciPy

## Examples

- [`examples/ds003416/`](examples/ds003416/) — MASiVar dataset: 5 subjects, 3 same-vendor sites (single-shell)
- [`examples/ds005664/`](examples/ds005664/) — SDSU Traveling Subjects: 9 subjects, cross-vendor GE vs Siemens (multi-shell)

## Documentation

- [Design Document](docs/design/DESIGN.md)
- [Pipeline Flowchart](docs/design/FLOWCHART.md)
- [Examples](examples/)

## Citation

If you use this tool, please cite:

**RISH harmonization method:**
- Mirzaalian H, et al. (2016). Inter-site and inter-scanner diffusion MRI data harmonization. *NeuroImage*, 135, 311-323.
- Cetin Karayumak S, et al. (2019). Retrospective harmonization of multi-site diffusion MRI data acquired with different acquisition parameters. *NeuroImage*, 184, 180-200.

**RISH-GLM joint modeling:**
- De Luca A, et al. (2025). RISH-GLM: Rotationally Invariant Spherical Harmonic General Linear Model for quantitative dMRI harmonization. *Magnetic Resonance in Medicine*.

**MRtrix3:**
- Tournier JD, et al. (2019). MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. *NeuroImage*, 202, 116137.

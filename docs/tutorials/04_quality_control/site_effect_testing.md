# Statistical Testing for Site Effects

This tutorial explains how to use `mrtrix-rish site-effect` to statistically test
whether site effects are present in your multi-site diffusion MRI data, and to
validate that harmonization successfully removed them.

## Background

In multi-site studies, scanner differences introduce systematic biases that
confound biological signal. The site effect test fits a voxel-wise General Linear
Model (GLM) with site as a factor and uses Freedman-Lane permutation testing to
produce robust, non-parametric p-values. The result tells you whether there is a
statistically significant difference between sites at each voxel.

**Typical workflow:**

1. Run the test on **pre-harmonization** data — expect a large site effect.
2. Harmonize with `mrtrix-rish harmonize`.
3. Run the test on **post-harmonization** data — expect the site effect to be gone.
4. Compare the two results with `compare_site_effects()`.

## Prerequisites

- MRtrix3 >= 3.0.4
- Python >= 3.9 with `mrtrix-rish` installed (`pip install -e ".[qc]"`)
- Scalar metric images (e.g. FA maps) for each subject, **registered to a common
  template space**
- A brain mask in the same template space

## Preparing the Site List CSV

The CLI expects a CSV file that maps each subject to a site and an image path.
The file **must** have `subject`, `site`, and an image path column (any of:
`image_path`, `image`, `path`, `fa_path`, or `fa`).

Optional covariate columns (e.g. `age`, `sex`) can be included and activated
with the `--covariates` flag.

**Example: `sites.csv`**

```csv
subject,site,image_path,age,sex
sub-01,SiteA,/data/template_space/sub-01_fa.mif,25,M
sub-02,SiteA,/data/template_space/sub-02_fa.mif,32,F
sub-03,SiteA,/data/template_space/sub-03_fa.mif,28,M
sub-04,SiteB,/data/template_space/sub-04_fa.mif,45,F
sub-05,SiteB,/data/template_space/sub-05_fa.mif,38,M
sub-06,SiteB,/data/template_space/sub-06_fa.mif,52,F
sub-07,SiteC,/data/template_space/sub-07_fa.mif,30,M
sub-08,SiteC,/data/template_space/sub-08_fa.mif,41,F
sub-09,SiteC,/data/template_space/sub-09_fa.mif,35,M
```

## Running the Test

### Basic usage (permutation test, recommended)

```bash
mrtrix-rish site-effect \
    --site-list sites.csv \
    --mask template_mask.mif \
    --output site_effect_pre/ \
    --n-permutations 5000 \
    --seed 42
```

### With covariates

Including covariates (e.g. age, sex) in the design matrix ensures that
biological variation is modelled as a nuisance regressor rather than being
conflated with the site effect.

```bash
mrtrix-rish site-effect \
    --site-list sites.csv \
    --mask template_mask.mif \
    --output site_effect_pre/ \
    --covariates age,sex \
    --n-permutations 5000 \
    --seed 42
```

Sex values are automatically coded as numeric (M/MALE/1 = 1.0, everything
else = 0.0). Continuous covariates are z-scored by default.

### Heteroscedastic test

If you suspect that different scanners produce different noise levels (common in
multi-vendor studies), use the `--heteroscedastic` flag. This switches from the
standard F-statistic to a Welch-type G-statistic that allows per-site variance
estimates.

```bash
mrtrix-rish site-effect \
    --site-list sites.csv \
    --mask template_mask.mif \
    --output site_effect_pre/ \
    --heteroscedastic \
    --n-permutations 5000
```

### All CLI options

| Flag | Default | Description |
|------|---------|-------------|
| `--site-list`, `-s` | (required) | CSV with subject, site, and image path columns |
| `--mask`, `-m` | (required) | Brain mask in template space |
| `--output`, `-o` | (required) | Output directory |
| `--test`, `-t` | `permutation` | Test type: `permutation` or `parametric` |
| `--n-permutations`, `-n` | `5000` | Number of permutations |
| `--covariates`, `-c` | None | Comma-separated covariate column names |
| `--alpha`, `-a` | `0.05` | Significance level |
| `--seed` | None | Random seed for reproducibility |
| `--heteroscedastic` | off | Use per-site variance (G-statistic) |

## Output Files

The command writes the following to the output directory:

```
site_effect_pre/
├── summary.json              # Key summary statistics
├── f_statistic.mif           # Voxel-wise F (or G) statistic map
├── p_permutation.mif         # Permutation p-value map
├── effect_size_eta_sq.mif    # Partial eta-squared effect size map
├── neg_log10_p.mif           # -log10(p) map for visualization
└── significant_mask.mif      # Binary mask of significant voxels (p < alpha)
```

### summary.json

```json
{
  "n_subjects": 150,
  "n_sites": 3,
  "n_voxels": 98432,
  "n_permutations": 5000,
  "unique_sites": ["SiteA", "SiteB", "SiteC"],
  "percent_significant_uncorrected": 45.2,
  "percent_significant_fdr": 38.7,
  "percent_significant_permutation": 42.1,
  "mean_effect_size": 0.12,
  "median_effect_size": 0.08,
  "fdr_threshold": 0.032
}
```

## Interpreting the Results

### Key metrics

- **percent_significant_permutation** — Percentage of voxels with a significant
  site effect (permutation p < alpha). This is the primary outcome measure.
- **mean_effect_size** — Mean partial eta-squared across voxels. Indicates the
  magnitude of the site effect (0 = none, 1 = total).
- **percent_significant_fdr** — Percentage after FDR correction (more
  conservative than uncorrected).

### Decision thresholds

| Metric | Pass | Warning | Fail |
|--------|------|---------|------|
| % significant voxels (perm.) | < 5% | 5-15% | > 15% |
| Mean effect size (eta-sq) | < 0.02 | 0.02-0.06 | > 0.06 |

### Viewing the maps

Use `mrview` (MRtrix3 viewer) to inspect the spatial distribution of site effects:

```bash
# Overlay -log10(p) on a template FA
mrview template_fa.mif \
    -overlay.load site_effect_pre/neg_log10_p.mif \
    -overlay.colourmap hot \
    -overlay.threshold_min 1.3

# View significant voxels
mrview template_fa.mif \
    -overlay.load site_effect_pre/significant_mask.mif \
    -overlay.colourmap red
```

The threshold of 1.3 on the `-log10(p)` map corresponds to p = 0.05.

## Full Example: Pre vs Post Harmonization

### Step 1: Test before harmonization

```bash
mrtrix-rish site-effect \
    --site-list sites_pre.csv \
    --mask template_mask.mif \
    --covariates age,sex \
    --output site_effect_pre/ \
    --n-permutations 5000 \
    --seed 42
```

### Step 2: Harmonize

```bash
mrtrix-rish create-template \
    --reference-list reference_subjects.txt \
    --output template/

mrtrix-rish harmonize \
    --target target_sh.mif \
    --template template/ \
    --output harmonized/
```

### Step 3: Test after harmonization

Regenerate the scalar maps (FA) from harmonized SH data, then run the same test:

```bash
mrtrix-rish site-effect \
    --site-list sites_post.csv \
    --mask template_mask.mif \
    --covariates age,sex \
    --output site_effect_post/ \
    --n-permutations 5000 \
    --seed 42
```

Using the same `--seed` ensures the same permutation sequence, making the
comparison more directly interpretable.

### Step 4: Compare results (Python API)

```python
from src.qc.site_effects import test_site_effect, compare_site_effects
import json

# Load previously saved results, or run test_site_effect() directly
# ...

comparison = compare_site_effects(pre_result, post_result, output_dir="comparison/")

print(json.dumps(comparison, indent=2))
```

**Example comparison output (`comparison/comparison.json`):**

```json
{
  "pre_harmonization": {
    "percent_significant": 42.1,
    "mean_effect_size": 0.12,
    "median_effect_size": 0.08
  },
  "post_harmonization": {
    "percent_significant": 3.1,
    "mean_effect_size": 0.015,
    "median_effect_size": 0.009
  },
  "reduction": {
    "percent_significant_reduction": 92.6,
    "effect_size_reduction": 87.5
  },
  "success_criteria": {
    "significant_voxels_below_5_percent": true,
    "effect_size_reduction_above_70_percent": true
  },
  "harmonization_successful": true
}
```

### Success criteria

The `compare_site_effects()` function applies two pass/fail checks:

1. **Post-harmonization significant voxels < 5%** — the residual site effect
   should be at or below the expected false-positive rate.
2. **Effect size reduction > 70%** — harmonization should substantially reduce
   the magnitude of the site effect.

Both must pass for `harmonization_successful` to be `true`.

## Statistical Details

### GLM formulation

The design matrix encodes site membership with dummy variables (treatment
coding, first site as reference) plus optional continuous/categorical covariates:

```
Y = intercept + site_B + site_C + ... + age + sex + error
```

The site effect hypothesis is an F-test across all site dummy variables,
testing H0: all site coefficients = 0.

### Permutation inference (Freedman-Lane)

Rather than assuming a parametric F-distribution, p-values are computed by
permuting subject labels and recomputing the F-statistic on each permutation.
The proportion of permutation statistics >= the observed statistic gives the
p-value.

The first "permutation" is always the identity (unpermuted data), so the
observed statistic is included in the null distribution.

### Multiple comparison correction

Two correction methods are applied:

- **FDR (Benjamini-Hochberg)** — controls the expected proportion of false
  positives among rejected voxels.
- **Max-statistic FWER** — available via the Python API (`max_statistic_correction()`),
  controls the family-wise error rate by comparing each voxel's statistic to the
  distribution of the maximum statistic across all voxels.

### Effect sizes

- **Partial eta-squared** — proportion of variance explained by the site factor
  after accounting for covariates. Range: 0 (no effect) to 1 (total effect).
  Computed as SS_between / (SS_between + SS_within).
- **Cohen's f** — derived from eta-squared: f = sqrt(eta_sq / (1 - eta_sq)).
  Convention: 0.10 = small, 0.25 = medium, 0.40 = large.

## Python API

For scripting and integration into larger pipelines, the full API is available:

```python
from src.qc.site_effects import (
    test_site_effect,
    compare_site_effects,
    fdr_correction,
    permutation_p_values,
    max_statistic_correction,
    compute_partial_eta_squared,
    compute_cohens_f,
    SiteEffectResult,
    Shuffler,
)
from src.qc.glm import (
    create_design_matrix,
    create_site_contrast,
    check_design,
    Hypothesis,
    TestFixedHomoscedastic,
    TestFixedHeteroscedastic,
)
```

### Custom hypothesis testing

You can define arbitrary contrasts beyond the built-in site effect test:

```python
import numpy as np
from src.qc.glm import Hypothesis, TestFixedHomoscedastic, create_design_matrix

# Example: test SiteB vs SiteC specifically (ignoring SiteA reference)
# Design columns: [intercept, site_B, site_C, age, sex]
contrast = np.array([[0, 1, -1, 0, 0]])  # t-test: site_B - site_C
hypothesis = Hypothesis(contrast, name="SiteB_vs_SiteC")

design, col_names = create_design_matrix(site_labels, covariates)
test = TestFixedHomoscedastic(data, design, [hypothesis])
result = test()[0]
# result.statistic contains the t-statistic per voxel
```

## Troubleshooting

**"High condition number" warning** — The design matrix is nearly
rank-deficient. This can happen when covariates are collinear with site
membership (e.g. all subjects at one site are the same sex). Consider removing
the problematic covariate or merging sites.

**All p-values are NaN** — scipy is not installed. Install it with
`pip install scipy` or use the permutation test (which does not require scipy
for the permutation p-values themselves, only for the parametric fallback).

**Very slow with many voxels** — Each permutation requires a full GLM solve
across all voxels. Reduce the mask to white matter only, or lower
`--n-permutations` to 1000 for a quick check before running the full 5000.

**No significant effect detected pre-harmonization** — If the pre-harmonization
test shows no site effect, harmonization may not be necessary. Verify that
images are correctly registered to template space and that the mask covers the
expected region.

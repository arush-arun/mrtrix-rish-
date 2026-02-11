# RISH Pipeline Scripts

Complete preprocessing + RISH harmonization pipeline for multi-site dMRI.

## Pipeline Steps

```
Raw DWI (BIDS)
    │
    ├── 1. Denoising         (dwidenoise)
    │
    ├── 2. Gibbs unringing   (mrdegibbs)
    │
    ├── 3. Distortion corr.  (dwifslpreproc - topup + eddy)
    │
    ├── 4. Bias correction   (dwibiascorrect ants)
    │
    ├── 5. Brain mask        (dwi2mask)
    │
    ├── 6. Response function (dwi2response tournier/dhollander)
    │
    ├── 7. FOD computation   (dwi2fod csd/msmt_csd)
    │
    ├── 8. RISH extraction   (θ₀, θ₂, θ₄, θ₆, θ₈)
    │
    ├── 9. Scale maps        (reference / target)
    │
    └── 10. Apply harmonization → Harmonized FOD
```

## Quick Start

1. **Edit configuration:**
   ```bash
   nano scripts/rish_config.conf
   ```
   
   Set:
   - `INPUT_PATH` - your BIDS dataset
   - `REFERENCE_SUBJECTS` - subjects from reference site
   - `TARGET_SUBJECTS` - subjects to harmonize

2. **Create subject list:**
   ```bash
   nano scripts/sublist.txt
   ```

3. **Run pipeline:**
   ```bash
   cd /home/uqahonne/uq/nif/mrtrix-rish
   ./scripts/rish_pipeline.sh
   ```

## Configuration Options

### Paths
```bash
INPUT_PATH=/path/to/bids          # BIDS dataset root
DERIVATIVES_PATH=derivatives/rish  # Output location
```

### BIDS Settings
```bash
SESSION=ses-1                      # Session to process
MAIN_PE_DIR=AP                     # Main phase encoding direction
BLIPPED_PE_DIR=PA                  # Reverse PE direction
```

### Processing Parameters
```bash
VOXEL_SIZE=1.25                    # Upsample voxel size (or "native")
LMAX=8                             # Maximum SH order
RESPONSE_ALGORITHM=auto            # auto, tournier, or dhollander
```

### RISH Harmonization
```bash
REFERENCE_SUBJECTS="sub-01 sub-02" # Reference site subjects
TARGET_SUBJECTS="sub-03 sub-04"    # Subjects to harmonize
SMOOTHING_FWHM=3.0                 # Scale map smoothing (mm)
SCALE_MIN=0.5                      # Minimum scale factor
SCALE_MAX=2.0                      # Maximum scale factor
```

### Processing Control
```bash
STEPS=all                          # all, or: preprocess,fod,rish,qc
ENABLE_CHECKPOINTS=true            # Resume capability
NTHREADS=4                         # Parallel threads
```

## Output Structure

```
derivatives/mrtrix-rish/
├── sub-XXX/
│   └── ses-1/
│       └── dwi/
│           ├── dwi_raw.mif
│           ├── dwi_denoised.mif
│           ├── dwi_unring.mif
│           ├── dwi_preproc.mif
│           ├── dwi_unbiased.mif
│           ├── dwi_final.mif
│           ├── mask.mif
│           ├── wm_fod.mif
│           ├── rish/
│           │   ├── rish_l0.mif
│           │   ├── rish_l2.mif
│           │   └── ...
│           └── harmonized/        # Target subjects only
│               ├── wm_fod_harmonized.mif
│               ├── scale_maps/
│               └── rish/
├── template/
│   ├── template_rish_l0.mif
│   └── ...
└── qc/
    └── harmonization_summary.csv
```

## Run Specific Steps

```bash
# Preprocessing only
STEPS=preprocess ./scripts/rish_pipeline.sh

# FOD only (requires preprocessing)
STEPS=fod ./scripts/rish_pipeline.sh

# RISH harmonization only (requires FOD)
STEPS=rish ./scripts/rish_pipeline.sh

# QC report only
STEPS=qc ./scripts/rish_pipeline.sh
```

## Resume Processing

With `ENABLE_CHECKPOINTS=true`, the pipeline creates checkpoint files in `./rish_checkpoints/`. 
Re-running the script will skip completed steps.

To force reprocessing:
```bash
rm -rf rish_checkpoints/
```

## References

- Cetin Karayumak et al. (2019). Retrospective harmonization of multi-site diffusion MRI data. NeuroImage.
- Tournier et al. (2019). MRtrix3: A fast, flexible and open software framework. NeuroImage.

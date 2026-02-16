#!/usr/bin/env bash
# =============================================================================
# End-to-end RISH harmonization on SDSU-TS (ds005664)
#
# Dataset: OpenNeuro ds005664 — SDSU Traveling Subjects
#   9 adults scanned on 2 cross-vendor scanners
#   Site SDSU: GE Discovery MR750 3T (ses-sdsu1)
#   Site CFMRI: Siemens Prisma 3T (ses-cfmri1)
#   Acquisition: 2-shell (b=1500/3000), 93 directions
#   Reference site: SDSU (GE)
#
# Requirements: MRtrix3 >= 3.0.8, mrtrix-rish (pip install -e ".[qc]")
#               AWS CLI (pip install awscli) for S3 download
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
NTHREADS="${NTHREADS:-4}"
N_PERMS="${N_PERMS:-1000}"

# Save system Python before modifying PATH
SYS_PYTHON="$(which python)"

# MRtrix3 conda environment
MRTRIX_ENV="${MRTRIX_ENV:-mrtrix3_env}"
MRTRIX_BIN="$HOME/anaconda3/envs/$MRTRIX_ENV/bin"
if [ -d "$MRTRIX_BIN" ]; then
    export PATH="$MRTRIX_BIN:$PATH"
    echo "Using MRtrix3 from: $MRTRIX_BIN ($(mrconvert --version 2>&1 | grep -oP 'mrconvert \K[0-9.]+' || echo 'unknown'))"
else
    echo "WARNING: $MRTRIX_BIN not found, using system MRtrix3"
fi
echo "Using Python for mrtrix-rish: $SYS_PYTHON"

# Directories
RAW_DIR="$SCRIPT_DIR/raw"
PROC_DIR="$SCRIPT_DIR/processing"
TEMPLATE_DIR="$SCRIPT_DIR/template"
TEMPLATE_COV_DIR="$SCRIPT_DIR/template_covariate"
HARM_DIR="$SCRIPT_DIR/harmonized"
HARM_COV_DIR="$SCRIPT_DIR/harmonized_covariate"
RESULTS_DIR="$SCRIPT_DIR/results"
PARTICIPANTS="$SCRIPT_DIR/participants.tsv"

# Reference scan for registration
REF_SUB="sub-ts001"
REF_SES="ses-sdsu1"

# ---- Subject manifest ----
# All 9 subjects (sub-ts006 does not exist in dataset)
SUBJECTS=("sub-ts001" "sub-ts002" "sub-ts003" "sub-ts004" "sub-ts005"
          "sub-ts007" "sub-ts008" "sub-ts009" "sub-ts010")

# Sessions per site
SES_SDSU="ses-sdsu1"
SES_CFMRI="ses-cfmri1"

# Acquisition label
ACQ="acq-2shell93dir"

# Subject demographics (for --subject-covariates)
declare -A SUBJ_AGE=(
    ["sub-ts001"]=38.7 ["sub-ts002"]=28.6 ["sub-ts003"]=22.4
    ["sub-ts004"]=39.9 ["sub-ts005"]=24.1 ["sub-ts007"]=39.9
    ["sub-ts008"]=23.7 ["sub-ts009"]=55.3 ["sub-ts010"]=25.3
)
declare -A SUBJ_SEX=(
    ["sub-ts001"]=F ["sub-ts002"]=F ["sub-ts003"]=F
    ["sub-ts004"]=M ["sub-ts005"]=F ["sub-ts007"]=M
    ["sub-ts008"]=F ["sub-ts009"]=F ["sub-ts010"]=M
)

# =============================================================================
# Step 1: Download raw DWI data from OpenNeuro
# =============================================================================
step1_download() {
    echo "============================================================"
    echo "Step 1: Downloading SDSU-TS data from OpenNeuro (ds005664)"
    echo "============================================================"

    if [ -d "$RAW_DIR" ] && [ "$(ls "$RAW_DIR"/sub-ts*/ses-*/dwi/*.nii.gz 2>/dev/null | wc -l)" -ge 18 ]; then
        echo "  [skip] Raw data already present ($(ls "$RAW_DIR"/sub-ts*/ses-*/dwi/*.nii.gz 2>/dev/null | wc -l) DWI files)"
        return
    fi

    mkdir -p "$RAW_DIR"

    echo "  Downloading 2-shell HCP-style DWI (acq-2shell93dir) for all subjects..."
    echo "  Sessions: ses-sdsu1 (GE) + ses-cfmri1 (Siemens)"
    echo "  Source: s3://openneuro.org/ds005664/"
    echo ""

    local S3_BASE="s3://openneuro.org/ds005664"

    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            local dwi_dir="$RAW_DIR/${sub}/${ses}/dwi"
            local prefix="${sub}_${ses}_${ACQ}_dir-AP_dwi"

            if [ -f "$dwi_dir/${prefix}.nii.gz" ]; then
                echo "  [skip] $sub $ses — already downloaded"
                continue
            fi

            echo "  Downloading $sub $ses ..."
            mkdir -p "$dwi_dir"

            for ext in nii.gz bval bvec json; do
                aws s3 cp --no-sign-request \
                    "${S3_BASE}/${sub}/${ses}/dwi/${prefix}.${ext}" \
                    "$dwi_dir/${prefix}.${ext}" \
                    --only-show-errors
            done
        done
    done

    echo ""
    n_files=$(ls "$RAW_DIR"/sub-ts*/ses-*/dwi/*.nii.gz 2>/dev/null | wc -l)
    echo "  Downloaded $n_files DWI files."
    echo "Download complete."
    echo ""
}

# =============================================================================
# Step 2: Convert to MIF, compute masks, response functions, and FODs
# =============================================================================
step2_process() {
    echo "============================================================"
    echo "Step 2: Converting, masking, and computing MSMT-CSD FODs"
    echo "============================================================"

    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            local dwi_dir="$RAW_DIR/${sub}/${ses}/dwi"
            local prefix="${sub}_${ses}_${ACQ}_dir-AP_dwi"
            local out="$PROC_DIR/${sub}/${ses}"

            if [ -f "$out/wm_fod.mif" ]; then
                echo "  [skip] $sub $ses — FOD exists"
                continue
            fi

            # Check that raw data exists
            if [ ! -f "$dwi_dir/${prefix}.nii.gz" ]; then
                echo "  [WARN] $sub $ses — no DWI file, skipping"
                continue
            fi

            echo "  Processing $sub $ses ..."
            mkdir -p "$out"

            # Convert NIfTI to MIF with gradient table
            mrconvert "$dwi_dir/${prefix}.nii.gz" "$out/dwi.mif" \
                -fslgrad "$dwi_dir/${prefix}.bvec" "$dwi_dir/${prefix}.bval" \
                -nthreads "$NTHREADS" -force -quiet

            # Brain mask
            dwi2mask "$out/dwi.mif" "$out/mask.mif" \
                -nthreads "$NTHREADS" -force -quiet

            # Multi-shell: 3-tissue response function estimation
            dwi2response dhollander "$out/dwi.mif" \
                "$out/response_wm.txt" "$out/response_gm.txt" "$out/response_csf.txt" \
                -mask "$out/mask.mif" \
                -nthreads "$NTHREADS" -force -quiet

            # MSMT-CSD FOD computation
            dwi2fod msmt_csd "$out/dwi.mif" \
                "$out/response_wm.txt" "$out/wm_fod.mif" \
                "$out/response_gm.txt" "$out/gm.mif" \
                "$out/response_csf.txt" "$out/csf.mif" \
                -mask "$out/mask.mif" -lmax 8,0,0 \
                -nthreads "$NTHREADS" -force -quiet

            echo "    Done: $out/wm_fod.mif"
        done
    done

    echo "FOD computation complete."
    echo ""
}

# =============================================================================
# Step 3: Register all FODs to a common reference
# =============================================================================
step3_register() {
    echo "============================================================"
    echo "Step 3: Registering FODs to common reference"
    echo "============================================================"

    local ref_dir="$PROC_DIR/${REF_SUB}/${REF_SES}"
    local ref_fod="$ref_dir/wm_fod.mif"
    local ref_mask="$ref_dir/mask.mif"

    echo "  Reference: $REF_SUB $REF_SES"

    # Copy reference as its own registered version
    cp -f "$ref_fod" "$ref_dir/fod_reg.mif" 2>/dev/null || true
    cp -f "$ref_mask" "$ref_dir/mask_reg.mif" 2>/dev/null || true

    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            local out="$PROC_DIR/${sub}/${ses}"

            # Skip if no FOD
            if [ ! -f "$out/wm_fod.mif" ]; then
                continue
            fi

            # Skip reference
            if [ "$sub" = "$REF_SUB" ] && [ "$ses" = "$REF_SES" ]; then
                continue
            fi

            if [ -f "$out/fod_reg.mif" ]; then
                echo "  [skip] $sub $ses — registered FOD exists"
                continue
            fi

            echo "  Registering $sub $ses ..."

            # Affine FOD registration
            mrregister "$out/wm_fod.mif" "$ref_fod" \
                -type affine \
                -affine "$out/affine.txt" \
                -mask1 "$out/mask.mif" -mask2 "$ref_mask" \
                -nthreads "$NTHREADS" -force -quiet

            # Transform FOD with reorientation
            mrtransform "$out/wm_fod.mif" "$out/fod_reg.mif" \
                -linear "$out/affine.txt" \
                -template "$ref_fod" \
                -reorient_fod yes \
                -nthreads "$NTHREADS" -force -quiet

            # Transform mask
            mrtransform "$out/mask.mif" "$out/mask_reg.mif" \
                -linear "$out/affine.txt" \
                -template "$ref_fod" \
                -interp nearest \
                -nthreads "$NTHREADS" -force -quiet

            echo "    Done."
        done
    done

    # Create group mask (intersection)
    echo "  Creating group mask (intersection) ..."
    local mask_list=""
    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            if [ -f "$PROC_DIR/${sub}/${ses}/mask_reg.mif" ]; then
                mask_list="$mask_list $PROC_DIR/${sub}/${ses}/mask_reg.mif"
            fi
        done
    done
    # shellcheck disable=SC2086
    mrmath $mask_list min "$PROC_DIR/group_mask.mif" -datatype bit -force -quiet

    echo "Registration complete."
    echo ""
}

# =============================================================================
# Step 4: Extract pre-harmonization RISH l0
# =============================================================================
step4_rish_pre() {
    echo "============================================================"
    echo "Step 4: Extracting pre-harmonization RISH l0 features"
    echo "============================================================"

    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            local out="$PROC_DIR/${sub}/${ses}"

            if [ ! -f "$out/fod_reg.mif" ]; then
                continue
            fi

            if [ -f "$out/rish_l0_pre.mif" ]; then
                echo "  [skip] $sub $ses"
                continue
            fi

            echo "  Extracting l0 for $sub $ses ..."
            mrconvert "$out/fod_reg.mif" -coord 3 0 "$out/rish_l0_pre.mif" -force -quiet
        done
    done

    echo "RISH l0 extraction complete."
    echo ""
}

# =============================================================================
# Step 5: Site-effect test BEFORE harmonization
# =============================================================================
step5_site_effect_pre() {
    echo "============================================================"
    echo "Step 5: Site-effect test (pre-harmonization, RISH l0)"
    echo "============================================================"

    local csv="$SCRIPT_DIR/sites_pre.csv"
    echo "subject,site,image_path" > "$csv"

    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            local rish="$PROC_DIR/${sub}/${ses}/rish_l0_pre.mif"
            if [ -f "$rish" ]; then
                local site
                if [ "$ses" = "$SES_SDSU" ]; then site="SDSU_GE"; else site="CFMRI_Siemens"; fi
                echo "${sub}_${ses},${site},${rish}" >> "$csv"
            fi
        done
    done

    mkdir -p "$RESULTS_DIR/pre"

    $SYS_PYTHON -m src.cli.main site-effect \
        --site-list "$csv" \
        --mask "$PROC_DIR/group_mask.mif" \
        --output "$RESULTS_DIR/pre" \
        --n-permutations "$N_PERMS" \
        --seed 42

    echo ""
}

# =============================================================================
# Step 6: RISH harmonization (baseline, no covariates)
# =============================================================================
step6_harmonize() {
    echo "============================================================"
    echo "Step 6: RISH harmonization (reference: SDSU/GE)"
    echo "============================================================"

    local group_mask="$PROC_DIR/group_mask.mif"

    # Create reference lists (SDSU subjects)
    local ref_list="$SCRIPT_DIR/sdsu_fods.txt"
    local mask_list="$SCRIPT_DIR/sdsu_masks.txt"
    > "$ref_list"
    > "$mask_list"

    for sub in "${SUBJECTS[@]}"; do
        local fod="$PROC_DIR/${sub}/${SES_SDSU}/fod_reg.mif"
        if [ -f "$fod" ]; then
            echo "$fod" >> "$ref_list"
            echo "$group_mask" >> "$mask_list"
        fi
    done

    # Create template
    echo "  Creating RISH template from SDSU (GE) subjects ..."
    mkdir -p "$TEMPLATE_DIR"

    $SYS_PYTHON -m src.cli.main create-template \
        --reference-list "$ref_list" \
        --mask-list "$mask_list" \
        --output "$TEMPLATE_DIR" \
        --lmax 8 \
        --nthreads "$NTHREADS"

    # Harmonize CFMRI (Siemens) subjects
    for sub in "${SUBJECTS[@]}"; do
        local fod="$PROC_DIR/${sub}/${SES_CFMRI}/fod_reg.mif"
        local out="$HARM_DIR/${sub}/${SES_CFMRI}"

        if [ ! -f "$fod" ]; then
            continue
        fi

        if [ -f "$out/sh_harmonized.mif" ]; then
            echo "  [skip] $sub $SES_CFMRI — harmonized"
            continue
        fi

        echo "  Harmonizing $sub $SES_CFMRI ..."
        mkdir -p "$out"

        $SYS_PYTHON -m src.cli.main harmonize \
            --target "$fod" \
            --mask "$group_mask" \
            --template "$TEMPLATE_DIR" \
            --output "$out" \
            --nthreads "$NTHREADS"
    done

    echo "Baseline harmonization complete."
    echo ""
}

# =============================================================================
# Step 7: Covariate-adjusted harmonization
# =============================================================================
step7_harmonize_covariate() {
    echo "============================================================"
    echo "Step 7: Covariate-adjusted RISH harmonization (age, sex)"
    echo "============================================================"

    local group_mask="$PROC_DIR/group_mask.mif"

    # Create reference lists (SDSU subjects)
    local ref_list="$SCRIPT_DIR/sdsu_fods.txt"
    local mask_list="$SCRIPT_DIR/sdsu_masks.txt"

    # Create covariate-adjusted template
    echo "  Creating covariate-adjusted template ..."
    mkdir -p "$TEMPLATE_COV_DIR"

    $SYS_PYTHON -m src.cli.main create-template \
        --reference-list "$ref_list" \
        --mask-list "$mask_list" \
        --output "$TEMPLATE_COV_DIR" \
        --lmax 8 \
        --nthreads "$NTHREADS" \
        --participants "$PARTICIPANTS" \
        --covariates age,sex \
        --subject-column participant_id

    # Harmonize CFMRI (Siemens) subjects with covariates
    for sub in "${SUBJECTS[@]}"; do
        local fod="$PROC_DIR/${sub}/${SES_CFMRI}/fod_reg.mif"
        local out="$HARM_COV_DIR/${sub}/${SES_CFMRI}"

        if [ ! -f "$fod" ]; then
            continue
        fi

        if [ -f "$out/sh_harmonized.mif" ]; then
            echo "  [skip] $sub $SES_CFMRI — covariate-harmonized"
            continue
        fi

        local age="${SUBJ_AGE[$sub]}"
        local sex="${SUBJ_SEX[$sub]}"

        echo "  Harmonizing $sub $SES_CFMRI (age=$age, sex=$sex) ..."
        mkdir -p "$out"

        $SYS_PYTHON -m src.cli.main harmonize \
            --target "$fod" \
            --mask "$group_mask" \
            --template "$TEMPLATE_COV_DIR" \
            --output "$out" \
            --nthreads "$NTHREADS" \
            --subject-covariates "age=${age},sex=${sex}"
    done

    echo "Covariate-adjusted harmonization complete."
    echo ""
}

# =============================================================================
# Step 8: Extract post-harmonization RISH l0
# =============================================================================
step8_rish_post() {
    echo "============================================================"
    echo "Step 8: Extracting post-harmonization RISH l0"
    echo "============================================================"

    # SDSU (reference): use pre-harmonization l0
    for sub in "${SUBJECTS[@]}"; do
        local src="$PROC_DIR/${sub}/${SES_SDSU}/rish_l0_pre.mif"
        if [ -f "$src" ]; then
            for suffix in post postcov; do
                local dst="$PROC_DIR/${sub}/${SES_SDSU}/rish_l0_${suffix}.mif"
                if [ ! -f "$dst" ]; then
                    cp "$src" "$dst"
                fi
            done
        fi
    done

    # CFMRI: extract l0 from baseline harmonized
    for sub in "${SUBJECTS[@]}"; do
        local harm_fod="$HARM_DIR/${sub}/${SES_CFMRI}/sh_harmonized.mif"
        local dst="$PROC_DIR/${sub}/${SES_CFMRI}/rish_l0_post.mif"

        if [ -f "$harm_fod" ] && [ ! -f "$dst" ]; then
            echo "  Extracting post (baseline) l0 for $sub $SES_CFMRI ..."
            mrconvert "$harm_fod" -coord 3 0 "$dst" -force -quiet
        fi
    done

    # CFMRI: extract l0 from covariate-adjusted harmonized
    for sub in "${SUBJECTS[@]}"; do
        local harm_fod="$HARM_COV_DIR/${sub}/${SES_CFMRI}/sh_harmonized.mif"
        local dst="$PROC_DIR/${sub}/${SES_CFMRI}/rish_l0_postcov.mif"

        if [ -f "$harm_fod" ] && [ ! -f "$dst" ]; then
            echo "  Extracting post (covariate) l0 for $sub $SES_CFMRI ..."
            mrconvert "$harm_fod" -coord 3 0 "$dst" -force -quiet
        fi
    done

    echo "Post-harmonization RISH l0 extraction complete."
    echo ""
}

# =============================================================================
# Step 9: Site-effect tests AFTER harmonization
# =============================================================================
step9_site_effect_post() {
    echo "============================================================"
    echo "Step 9: Site-effect tests (post-harmonization)"
    echo "============================================================"

    local group_mask="$PROC_DIR/group_mask.mif"

    # --- Baseline post ---
    local csv_post="$SCRIPT_DIR/sites_post.csv"
    echo "subject,site,image_path" > "$csv_post"

    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            local rish="$PROC_DIR/${sub}/${ses}/rish_l0_post.mif"
            if [ -f "$rish" ]; then
                local site
                if [ "$ses" = "$SES_SDSU" ]; then site="SDSU_GE"; else site="CFMRI_Siemens"; fi
                echo "${sub}_${ses},${site},${rish}" >> "$csv_post"
            fi
        done
    done

    mkdir -p "$RESULTS_DIR/post"
    echo "  Running site-effect test (baseline post-harmonization) ..."

    $SYS_PYTHON -m src.cli.main site-effect \
        --site-list "$csv_post" \
        --mask "$group_mask" \
        --output "$RESULTS_DIR/post" \
        --n-permutations "$N_PERMS" \
        --seed 42

    # --- Covariate post ---
    local csv_postcov="$SCRIPT_DIR/sites_postcov.csv"
    echo "subject,site,image_path" > "$csv_postcov"

    for sub in "${SUBJECTS[@]}"; do
        for ses in "$SES_SDSU" "$SES_CFMRI"; do
            local rish="$PROC_DIR/${sub}/${ses}/rish_l0_postcov.mif"
            if [ -f "$rish" ]; then
                local site
                if [ "$ses" = "$SES_SDSU" ]; then site="SDSU_GE"; else site="CFMRI_Siemens"; fi
                echo "${sub}_${ses},${site},${rish}" >> "$csv_postcov"
            fi
        done
    done

    mkdir -p "$RESULTS_DIR/postcov"
    echo ""
    echo "  Running site-effect test (covariate-adjusted post-harmonization) ..."

    $SYS_PYTHON -m src.cli.main site-effect \
        --site-list "$csv_postcov" \
        --mask "$group_mask" \
        --output "$RESULTS_DIR/postcov" \
        --n-permutations "$N_PERMS" \
        --seed 42

    echo ""
}

# =============================================================================
# Step 10: Compare pre vs post-baseline vs post-covariate
# =============================================================================
step10_compare() {
    echo "============================================================"
    echo "Step 10: Comparing pre vs post harmonization"
    echo "============================================================"

    mkdir -p "$RESULTS_DIR/comparison"

    $SYS_PYTHON -c "
import json
from pathlib import Path

pre = json.loads(Path('$RESULTS_DIR/pre/summary.json').read_text())
post = json.loads(Path('$RESULTS_DIR/post/summary.json').read_text())
postcov = json.loads(Path('$RESULTS_DIR/postcov/summary.json').read_text())

print('SDSU-TS (ds005664): GE Discovery MR750 vs Siemens Prisma')
print('=' * 62)
print()
print('                         % Sig Voxels   Mean Effect Size')
print('                         -----------   ----------------')
print(f'  Pre-harmonization:      {pre[\"percent_significant_permutation\"]:>8.1f}%      {pre[\"mean_effect_size\"]:.4f}')
print(f'  Post (baseline):        {post[\"percent_significant_permutation\"]:>8.1f}%      {post[\"mean_effect_size\"]:.4f}')
print(f'  Post (covariate-adj):   {postcov[\"percent_significant_permutation\"]:>8.1f}%      {postcov[\"mean_effect_size\"]:.4f}')
print()

sig_pre = pre['percent_significant_permutation']
sig_post = post['percent_significant_permutation']
sig_postcov = postcov['percent_significant_permutation']
es_pre = pre['mean_effect_size']
es_post = post['mean_effect_size']
es_postcov = postcov['mean_effect_size']

sig_reduction = (sig_pre - sig_post) / (sig_pre + 1e-10) * 100
es_reduction = (es_pre - es_post) / (es_pre + 1e-10) * 100

print(f'  Baseline reduction:')
print(f'    Significant voxels: {sig_reduction:.1f}% reduction')
print(f'    Effect size: {es_reduction:.1f}% reduction')
print()

comparison = {
    'dataset': 'SDSU-TS (ds005664)',
    'sites': {'reference': 'SDSU (GE Discovery MR750)', 'target': 'CFMRI (Siemens Prisma)'},
    'pre_harmonization': {
        'percent_significant': sig_pre,
        'mean_effect_size': es_pre,
    },
    'post_baseline': {
        'percent_significant': sig_post,
        'mean_effect_size': es_post,
    },
    'post_covariate_adjusted': {
        'percent_significant': sig_postcov,
        'mean_effect_size': es_postcov,
    },
    'reduction': {
        'percent_significant_reduction': round(sig_reduction, 1),
        'effect_size_reduction': round(es_reduction, 1),
    },
    'success_criteria': {
        'significant_voxels_below_5_percent': sig_post < 5.0,
        'effect_size_reduction_above_70_percent': es_reduction > 70,
    },
    'harmonization_successful': sig_post < 5.0 and es_reduction > 70,
}

Path('$RESULTS_DIR/comparison/comparison.json').write_text(
    json.dumps(comparison, indent=2)
)

print(json.dumps(comparison, indent=2))
print()
if comparison['harmonization_successful']:
    print('RESULT: Cross-vendor harmonization SUCCESSFUL')
else:
    print('RESULT: Harmonization needs review')
    if not comparison['success_criteria']['significant_voxels_below_5_percent']:
        print(f'  - Significant voxels still at {sig_post:.1f}% (target: <5%)')
    if not comparison['success_criteria']['effect_size_reduction_above_70_percent']:
        print(f'  - Effect size reduction only {es_reduction:.1f}% (target: >70%)')
"
}

# =============================================================================
# Main
# =============================================================================
main() {
    echo "================================================================"
    echo "mrtrix-rish end-to-end test: SDSU-TS ds005664"
    echo "  9 subjects x 2 cross-vendor scanners (GE vs Siemens)"
    echo "  Acquisition: 2-shell (b=1500/3000), 93 directions"
    echo "  Reference site: SDSU (GE Discovery MR750)"
    echo "  Threads: $NTHREADS | Permutations: $N_PERMS"
    echo "================================================================"
    echo ""

    cd "$PROJECT_DIR"

    step1_download
    step2_process
    step3_register
    step4_rish_pre
    step5_site_effect_pre
    step6_harmonize
    step7_harmonize_covariate
    step8_rish_post
    step9_site_effect_post
    step10_compare

    echo ""
    echo "================================================================"
    echo "Pipeline complete. Results in: $RESULTS_DIR/"
    echo "================================================================"
}

show_help() {
    cat <<'HELP'
End-to-end RISH harmonization test on SDSU-TS (ds005664)

Usage: run_harmonization_test.sh [STEP] [OPTIONS]

Run all steps:
  run_harmonization_test.sh          Run the full pipeline (steps 1-10)
  run_harmonization_test.sh all      Same as above

Run individual steps:
  1  | download     Download DWI from OpenNeuro S3 via AWS CLI (~2-3 GB)
  2  | process      Convert to MIF, compute brain masks and MSMT-CSD FODs
  3  | register     Register all FODs to reference (sub-ts001/ses-sdsu1),
                    compute group mask
  4  | rish-pre     Extract pre-harmonization RISH l0 from registered FODs
  5  | site-pre     Run GLM site-effect test on pre-harmonization RISH l0
  6  | harmonize    Create RISH template from SDSU (GE), harmonize CFMRI
                    (Siemens) to match
  7  | harm-cov     Covariate-adjusted harmonization (age, sex)
  8  | rish-post    Extract post-harmonization RISH l0
  9  | site-post    Run site-effect tests on post-harmonization data
  10 | compare      Compare pre vs post results, report success/failure

Environment variables:
  NTHREADS     Number of threads for MRtrix3 commands (default: 4)
  N_PERMS      Number of permutations for site-effect test (default: 1000)
  MRTRIX_ENV   Conda environment with MRtrix3 (default: mrtrix3_env)

Requirements:
  - MRtrix3 >= 3.0.8 (via conda env)
  - mrtrix-rish with QC extras: pip install -e ".[qc]"
  - AWS CLI: pip install awscli (for S3 download)
  - ~10-15 GB disk space

Dataset:
  OpenNeuro ds005664 (SDSU Traveling Subjects)
  9 adults x 2 cross-vendor scanners
  GE Discovery MR750 (SDSU) vs Siemens Prisma (CFMRI)
  2-shell: b=1500/3000, 93 directions

Examples:
  bash run_harmonization_test.sh                # full pipeline
  bash run_harmonization_test.sh download       # just download data
  bash run_harmonization_test.sh 6              # just run harmonization
  NTHREADS=8 bash run_harmonization_test.sh all # full pipeline, 8 threads
HELP
}

# Allow running individual steps
if [ $# -gt 0 ]; then
    case "$1" in
        -h|--help|help) show_help ;;
        1|download)     step1_download ;;
        2|process)      step2_process ;;
        3|register)     step3_register ;;
        4|rish-pre)     step4_rish_pre ;;
        5|site-pre)     step5_site_effect_pre ;;
        6|harmonize)    step6_harmonize ;;
        7|harm-cov)     step7_harmonize_covariate ;;
        8|rish-post)    step8_rish_post ;;
        9|site-post)    step9_site_effect_post ;;
        10|compare)     step10_compare ;;
        all)            main ;;
        *)              echo "Unknown step: $1"; echo ""; show_help; exit 1 ;;
    esac
else
    main
fi

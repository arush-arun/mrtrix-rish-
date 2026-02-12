#!/usr/bin/env bash
# =============================================================================
# End-to-end RISH harmonization + site-effect test on MASiVar (ds003416)
#
# Dataset: OpenNeuro ds003416 — MASiVar Cohort II
#   5 adults scanned on 3 scanners (s1A, s1B, s3)
#   Acquisition: b=1000, 96 directions
#   Reference site: s1B
#
# Requirements: MRtrix3 >= 3.0.8, mrtrix-rish (pip install -e ".[qc]")
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BASE_URL="https://s3.amazonaws.com/openneuro.org/ds003416"
NTHREADS="${NTHREADS:-4}"
N_PERMS="${N_PERMS:-1000}"

# Save system Python (for mrtrix-rish CLI) before modifying PATH
SYS_PYTHON="$(which python)"

# MRtrix3 conda environment — prepend to PATH so both bash and Python subprocess
# calls pick up the correct MRtrix3 binaries (set MRTRIX_ENV to override)
MRTRIX_ENV="${MRTRIX_ENV:-mrtrix3_env}"
MRTRIX_BIN="$HOME/anaconda3/envs/$MRTRIX_ENV/bin"
if [ -d "$MRTRIX_BIN" ]; then
    export PATH="$MRTRIX_BIN:$PATH"
    echo "Using MRtrix3 from: $MRTRIX_BIN ($(mrconvert --version 2>&1 | head -1))"
else
    echo "WARNING: $MRTRIX_BIN not found, using system MRtrix3"
fi
echo "Using Python for mrtrix-rish: $SYS_PYTHON"

# Directories
RAW_DIR="$SCRIPT_DIR/raw"
PROC_DIR="$SCRIPT_DIR/processing"
TEMPLATE_DIR="$SCRIPT_DIR/template"
HARM_DIR="$SCRIPT_DIR/harmonized"
RESULTS_DIR="$SCRIPT_DIR/results"

# Reference scan for registration
REF_SUB="sub-cIIs00"
REF_SES="ses-s1Bx1"

# ---- File manifest (subject session acq_label run_number) ----
# s1A files
declare -a S1A_FILES=(
    "sub-cIIs00 ses-s1Ax1 b1000n96r25x25x25peAPP 110"
    "sub-cIIs01 ses-s1Ax1 b1000n96r19x19x22peAPP 109"
    "sub-cIIs02 ses-s1Ax1 b1000n96r19x19x22peAPP 110"
    "sub-cIIs03 ses-s1Ax1 b1000n96r19x19x22peAPP 109"
    "sub-cIIs04 ses-s1Ax1 b1000n96r19x19x22peAPP 109"
)

# s1B files
declare -a S1B_FILES=(
    "sub-cIIs00 ses-s1Bx1 b1000n96r25x25x25peAPP 109"
    "sub-cIIs01 ses-s1Bx1 b1000n96r25x25x25peAPP 109"
    "sub-cIIs02 ses-s1Bx1 b1000n96r25x25x25peAPP 109"
    "sub-cIIs03 ses-s1Bx1 b1000n96r25x25x25peAPP 109"
    "sub-cIIs04 ses-s1Bx1 b1000n96r25x25x25peAPP 109"
)

# s3 files
declare -a S3_FILES=(
    "sub-cIIs00 ses-s3x1 b1000n96r25x25x25peAPP 110"
    "sub-cIIs01 ses-s3x1 b1000n96r25x25x25peAPP 110"
    "sub-cIIs02 ses-s3x1 b1000n96r25x25x25peAPP 109"
    "sub-cIIs03 ses-s3x1 b1000n96r25x25x25peAPP 110"
    "sub-cIIs04 ses-s3x1 b1000n96r25x25x25peAPP 110"
)

# Combine all files
declare -a ALL_FILES=("${S1A_FILES[@]}" "${S1B_FILES[@]}" "${S3_FILES[@]}")

# Map session prefix to site label
site_from_session() {
    case "$1" in
        ses-s1Ax*) echo "s1A" ;;
        ses-s1Bx*) echo "s1B" ;;
        ses-s3x*)  echo "s3"  ;;
        *)         echo "unknown" ;;
    esac
}

# =============================================================================
# Step 1: Download raw DWI data from OpenNeuro S3
# =============================================================================
step1_download() {
    echo "============================================================"
    echo "Step 1: Downloading raw DWI data from OpenNeuro"
    echo "============================================================"

    for entry in "${ALL_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local prefix="${sub}_${ses}_acq-${acq}_run-${run}_dwi"
        local remote_dir="${sub}/${ses}/dwi"
        local local_dir="$RAW_DIR/${sub}/${ses}"
        mkdir -p "$local_dir"

        for ext in nii.gz bval bvec; do
            local url="${BASE_URL}/${remote_dir}/${prefix}.${ext}"
            local dest="$local_dir/${prefix}.${ext}"
            if [ -f "$dest" ]; then
                echo "  [skip] $dest (exists)"
            else
                echo "  [download] ${prefix}.${ext}"
                curl -sL "$url" -o "$dest"
            fi
        done
    done

    echo "Download complete."
    echo ""
}

# =============================================================================
# Step 2 & 3: Convert to MIF, compute masks, response functions, and FODs
# =============================================================================
step2_3_process() {
    echo "============================================================"
    echo "Steps 2-3: Converting, masking, and computing FODs"
    echo "============================================================"

    for entry in "${ALL_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local prefix="${sub}_${ses}_acq-${acq}_run-${run}_dwi"
        local raw="$RAW_DIR/${sub}/${ses}"
        local out="$PROC_DIR/${sub}/${ses}"
        mkdir -p "$out"

        if [ -f "$out/fod.mif" ]; then
            echo "  [skip] $sub $ses — FOD exists"
            continue
        fi

        echo "  Processing $sub $ses ..."

        # Convert NIfTI to MIF with gradient table
        mrconvert "$raw/${prefix}.nii.gz" "$out/dwi.mif" \
            -fslgrad "$raw/${prefix}.bvec" "$raw/${prefix}.bval" \
            -nthreads "$NTHREADS" -force -quiet

        # Brain mask
        dwi2mask "$out/dwi.mif" "$out/mask.mif" \
            -nthreads "$NTHREADS" -force -quiet

        # Response function estimation
        dwi2response tournier "$out/dwi.mif" "$out/response.txt" \
            -mask "$out/mask.mif" \
            -nthreads "$NTHREADS" -force -quiet

        # CSD FOD computation
        dwi2fod csd "$out/dwi.mif" "$out/response.txt" "$out/fod.mif" \
            -mask "$out/mask.mif" -lmax 8 \
            -nthreads "$NTHREADS" -force -quiet

        echo "    Done: $out/fod.mif"
    done

    echo "FOD computation complete."
    echo ""
}

# =============================================================================
# Step 4: Register all FODs to a common reference
# =============================================================================
step4_register() {
    echo "============================================================"
    echo "Step 4: Registering FODs to common reference"
    echo "============================================================"

    local ref_dir="$PROC_DIR/${REF_SUB}/${REF_SES}"
    local ref_fod="$ref_dir/fod.mif"
    local ref_mask="$ref_dir/mask.mif"

    echo "  Reference: $REF_SUB $REF_SES"

    # Copy reference as its own registered version
    local ref_reg="$PROC_DIR/${REF_SUB}/${REF_SES}"
    cp -f "$ref_fod" "$ref_reg/fod_reg.mif" 2>/dev/null || true
    cp -f "$ref_mask" "$ref_reg/mask_reg.mif" 2>/dev/null || true
    cp -f "$ref_dir/dwi.mif" "$ref_reg/dwi_reg.mif" 2>/dev/null || true

    for entry in "${ALL_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local out="$PROC_DIR/${sub}/${ses}"

        # Skip reference (already in target space)
        if [ "$sub" = "$REF_SUB" ] && [ "$ses" = "$REF_SES" ]; then
            continue
        fi

        if [ -f "$out/fod_reg.mif" ]; then
            echo "  [skip] $sub $ses — registered FOD exists"
            continue
        fi

        echo "  Registering $sub $ses ..."

        # Affine FOD registration
        mrregister "$out/fod.mif" "$ref_fod" \
            -type affine \
            -affine "$out/affine.txt" \
            -mask1 "$out/mask.mif" -mask2 "$ref_mask" \
            -nthreads "$NTHREADS" -force -quiet

        # Transform FOD (with reorientation for SH data)
        mrtransform "$out/fod.mif" "$out/fod_reg.mif" \
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

        # Transform DWI (for FA computation later)
        mrtransform "$out/dwi.mif" "$out/dwi_reg.mif" \
            -linear "$out/affine.txt" \
            -template "$ref_fod" \
            -nthreads "$NTHREADS" -force -quiet

        echo "    Done."
    done

    # Create group mask (intersection of all registered masks)
    echo "  Creating group mask (intersection) ..."
    local mask_list=""
    for entry in "${ALL_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        mask_list="$mask_list $PROC_DIR/${sub}/${ses}/mask_reg.mif"
    done
    # shellcheck disable=SC2086
    mrmath $mask_list min "$PROC_DIR/group_mask.mif" -datatype bit -force -quiet

    echo "Registration complete."
    echo ""
}

# =============================================================================
# Step 5: Compute pre-harmonization RISH features (l0 = isotropic FOD component)
# =============================================================================
step5_compute_rish_pre() {
    echo "============================================================"
    echo "Step 5: Extracting pre-harmonization RISH l0 features"
    echo "============================================================"

    for entry in "${ALL_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local out="$PROC_DIR/${sub}/${ses}"

        if [ -f "$out/rish_l0_pre.mif" ]; then
            echo "  [skip] $sub $ses — RISH l0 exists"
            continue
        fi

        echo "  Extracting RISH l0 for $sub $ses ..."

        # l0 = DC component of FOD (volume 0 of SH coefficients)
        mrconvert "$out/fod_reg.mif" -coord 3 0 "$out/rish_l0_pre.mif" -force -quiet

        echo "    Done."
    done

    echo "RISH l0 extraction complete."
    echo ""
}

# =============================================================================
# Step 6: Run site-effect test BEFORE harmonization
# =============================================================================
step6_site_effect_pre() {
    echo "============================================================"
    echo "Step 6: Site-effect test (pre-harmonization, RISH l0)"
    echo "============================================================"

    local csv="$SCRIPT_DIR/sites_pre.csv"
    echo "subject,site,image_path" > "$csv"

    for entry in "${ALL_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local site
        site=$(site_from_session "$ses")
        local rish="$PROC_DIR/${sub}/${ses}/rish_l0_pre.mif"
        echo "${sub}_${ses},${site},${rish}" >> "$csv"
    done

    mkdir -p "$RESULTS_DIR/pre"

    $SYS_PYTHON -m src.cli.main site-effect \
        --site-list "$csv" \
        --mask "$PROC_DIR/group_mask.mif" \
        --output "$RESULTS_DIR/pre" \
        --n-permutations "$N_PERMS" \
        --seed 42

    echo ""
    echo "Pre-harmonization results:"
    cat "$RESULTS_DIR/pre/summary.json"
    echo ""
}

# =============================================================================
# Step 7: RISH harmonization
# =============================================================================
step7_harmonize() {
    echo "============================================================"
    echo "Step 7: RISH harmonization (reference: s1B)"
    echo "============================================================"

    local group_mask="$PROC_DIR/group_mask.mif"

    # Create reference list (s1B FODs)
    local ref_list="$SCRIPT_DIR/s1b_fods.txt"
    local mask_list="$SCRIPT_DIR/s1b_masks.txt"
    > "$ref_list"
    > "$mask_list"

    for entry in "${S1B_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        echo "$PROC_DIR/${sub}/${ses}/fod_reg.mif" >> "$ref_list"
        echo "$group_mask" >> "$mask_list"
    done

    # Create RISH template from reference site
    echo "  Creating RISH template from s1B subjects ..."
    mkdir -p "$TEMPLATE_DIR"

    $SYS_PYTHON -m src.cli.main create-template \
        --reference-list "$ref_list" \
        --mask-list "$mask_list" \
        --output "$TEMPLATE_DIR" \
        --lmax 8 \
        --nthreads "$NTHREADS"

    # Harmonize target site subjects (s1A and s3)
    for entry in "${S1A_FILES[@]}" "${S3_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local out="$HARM_DIR/${sub}/${ses}"

        if [ -f "$out/sh_harmonized.mif" ]; then
            echo "  [skip] $sub $ses — harmonized FOD exists"
            continue
        fi

        echo "  Harmonizing $sub $ses ..."
        mkdir -p "$out"

        $SYS_PYTHON -m src.cli.main harmonize \
            --target "$PROC_DIR/${sub}/${ses}/fod_reg.mif" \
            --mask "$group_mask" \
            --template "$TEMPLATE_DIR" \
            --output "$out" \
            --nthreads "$NTHREADS"
    done

    echo "Harmonization complete."
    echo ""
}

# =============================================================================
# Step 8: Extract post-harmonization RISH l0 features
# =============================================================================
step8_compute_rish_post() {
    echo "============================================================"
    echo "Step 8: Extracting post-harmonization RISH l0 features"
    echo "============================================================"

    # For reference site (s1B): RISH l0 is unchanged (they define the template)
    for entry in "${S1B_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local src="$PROC_DIR/${sub}/${ses}/rish_l0_pre.mif"
        local dst="$PROC_DIR/${sub}/${ses}/rish_l0_post.mif"
        if [ ! -f "$dst" ]; then
            cp "$src" "$dst"
        fi
    done

    # For target sites: extract l0 from harmonized FOD
    for entry in "${S1A_FILES[@]}" "${S3_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local harm_dir="$HARM_DIR/${sub}/${ses}"
        local proc_dir="$PROC_DIR/${sub}/${ses}"

        if [ -f "$proc_dir/rish_l0_post.mif" ]; then
            echo "  [skip] $sub $ses — post RISH l0 exists"
            continue
        fi

        echo "  Extracting post-harmonization RISH l0 for $sub $ses ..."

        # l0 = DC component (volume 0) from harmonized FOD
        mrconvert "$harm_dir/sh_harmonized.mif" -coord 3 0 \
            "$proc_dir/rish_l0_post.mif" -force -quiet

        echo "    Done."
    done

    echo "Post-harmonization RISH l0 extraction complete."
    echo ""
}

# =============================================================================
# Step 9: Run site-effect test AFTER harmonization
# =============================================================================
step9_site_effect_post() {
    echo "============================================================"
    echo "Step 9: Site-effect test (post-harmonization, RISH l0)"
    echo "============================================================"

    local csv="$SCRIPT_DIR/sites_post.csv"
    echo "subject,site,image_path" > "$csv"

    for entry in "${ALL_FILES[@]}"; do
        read -r sub ses acq run <<< "$entry"
        local site
        site=$(site_from_session "$ses")
        local rish="$PROC_DIR/${sub}/${ses}/rish_l0_post.mif"
        echo "${sub}_${ses},${site},${rish}" >> "$csv"
    done

    mkdir -p "$RESULTS_DIR/post"

    $SYS_PYTHON -m src.cli.main site-effect \
        --site-list "$csv" \
        --mask "$PROC_DIR/group_mask.mif" \
        --output "$RESULTS_DIR/post" \
        --n-permutations "$N_PERMS" \
        --seed 42

    echo ""
    echo "Post-harmonization results:"
    cat "$RESULTS_DIR/post/summary.json"
    echo ""
}

# =============================================================================
# Step 10: Compare pre vs post
# =============================================================================
step10_compare() {
    echo "============================================================"
    echo "Step 10: Comparing pre vs post harmonization"
    echo "============================================================"

    mkdir -p "$RESULTS_DIR/comparison"

    python -c "
import sys
sys.path.insert(0, '$PROJECT_DIR')

import json
from pathlib import Path

pre = json.loads(Path('$RESULTS_DIR/pre/summary.json').read_text())
post = json.loads(Path('$RESULTS_DIR/post/summary.json').read_text())

sig_pre = pre['percent_significant_permutation']
sig_post = post['percent_significant_permutation']
es_pre = pre['mean_effect_size']
es_post = post['mean_effect_size']

sig_reduction = (sig_pre - sig_post) / (sig_pre + 1e-10) * 100
es_reduction = (es_pre - es_post) / (es_pre + 1e-10) * 100

comparison = {
    'pre_harmonization': {
        'percent_significant': sig_pre,
        'mean_effect_size': es_pre,
    },
    'post_harmonization': {
        'percent_significant': sig_post,
        'mean_effect_size': es_post,
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
    print('RESULT: Harmonization SUCCESSFUL')
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
    echo "mrtrix-rish end-to-end test: MASiVar ds003416 Cohort II"
    echo "  5 subjects x 3 scanners (s1A, s1B, s3)"
    echo "  Reference site: s1B"
    echo "  Threads: $NTHREADS | Permutations: $N_PERMS"
    echo "================================================================"
    echo ""

    cd "$PROJECT_DIR"

    step1_download
    step2_3_process
    step4_register
    step5_compute_rish_pre
    step6_site_effect_pre
    step7_harmonize
    step8_compute_rish_post
    step9_site_effect_post
    step10_compare

    echo ""
    echo "================================================================"
    echo "Pipeline complete. Results in: $RESULTS_DIR/"
    echo "================================================================"
}

show_help() {
    cat <<'HELP'
End-to-end RISH harmonization + site-effect test on MASiVar (ds003416)

Usage: run_site_effect_test.sh [STEP] [OPTIONS]

Run all steps:
  run_site_effect_test.sh          Run the full pipeline (steps 1-10)
  run_site_effect_test.sh all      Same as above

Run individual steps:
  1  | download     Download raw DWI from OpenNeuro S3 (~489 MB)
  2  | process      Convert to MIF, compute brain masks and CSD FODs (lmax=8)
  3  | fod          Alias for step 2
  4  | register     Register all FODs to reference (sub-cIIs00/ses-s1Bx1),
                    compute group mask
  5  | rish-pre     Extract pre-harmonization RISH l0 from registered FODs
  6  | site-pre     Run GLM site-effect test on pre-harmonization RISH l0
  7  | harmonize    Create RISH template from reference site (s1B),
                    harmonize target sites (s1A, s3) to match
  8  | rish-post    Extract post-harmonization RISH l0 from harmonized FODs
  9  | site-post    Run GLM site-effect test on post-harmonization RISH l0
  10 | compare      Compare pre vs post results, report success/failure

Environment variables:
  NTHREADS     Number of threads for MRtrix3 commands (default: 4)
  N_PERMS      Number of permutations for site-effect test (default: 1000)
  MRTRIX_ENV   Conda environment with MRtrix3 (default: mrtrix3_env)

Requirements:
  - MRtrix3 >= 3.0.8 (via conda env)
  - mrtrix-rish with QC extras: pip install -e ".[qc]"
  - ~2 GB disk space (489 MB raw + processing intermediates)

Dataset:
  OpenNeuro ds003416 (MASiVar) — Cohort II
  5 adults x 3 scanners (s1A, s1B, s3), b=1000, 96 directions
  Reference site: s1B

Examples:
  bash run_site_effect_test.sh                # full pipeline
  bash run_site_effect_test.sh download       # just download data
  bash run_site_effect_test.sh 7              # just run harmonization
  NTHREADS=8 bash run_site_effect_test.sh all # full pipeline, 8 threads
HELP
}

# Allow running individual steps
if [ $# -gt 0 ]; then
    case "$1" in
        -h|--help|help) show_help ;;
        1|download)     step1_download ;;
        2|process)      step2_3_process ;;
        3|fod)          step2_3_process ;;
        4|register)     step4_register ;;
        5|rish-pre)     step5_compute_rish_pre ;;
        6|site-pre)     step6_site_effect_pre ;;
        7|harmonize)    step7_harmonize ;;
        8|rish-post)    step8_compute_rish_post ;;
        9|site-post)    step9_site_effect_post ;;
        10|compare)     step10_compare ;;
        all)            main ;;
        *)              echo "Unknown step: $1"; echo ""; show_help; exit 1 ;;
    esac
else
    main
fi

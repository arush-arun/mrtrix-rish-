#!/usr/bin/env bash
# =============================================================================
# Phase 2 Test: Covariate-adjusted RISH harmonization on MASiVar (ds003416)
#
# Tests the covariate regression pipeline using age + sex as covariates.
# Compares results against baseline (no covariate adjustment).
#
# Prerequisite: steps 1-4 of run_site_effect_test.sh already completed
#   (download, process FODs, register to common space)
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
NTHREADS="${NTHREADS:-4}"
N_PERMS="${N_PERMS:-1000}"

SYS_PYTHON="$(which python)"

MRTRIX_ENV="${MRTRIX_ENV:-mrtrix3_env}"
MRTRIX_BIN="$HOME/anaconda3/envs/$MRTRIX_ENV/bin"
if [ -d "$MRTRIX_BIN" ]; then
    export PATH="$MRTRIX_BIN:$PATH"
    echo "Using MRtrix3 from: $MRTRIX_BIN"
fi
echo "Using Python: $SYS_PYTHON"

PROC_DIR="$SCRIPT_DIR/processing"
TEMPLATE_COV_DIR="$SCRIPT_DIR/template_covariate"
HARM_COV_DIR="$SCRIPT_DIR/harmonized_covariate"
RESULTS_DIR="$SCRIPT_DIR/results"
PARTICIPANTS="$SCRIPT_DIR/ref_participants.tsv"

# Subjects and sessions per site
declare -a S1B_SUBS=("sub-cIIs00" "sub-cIIs01" "sub-cIIs02" "sub-cIIs03" "sub-cIIs04")
declare -a S1B_SES=("ses-s1Bx1" "ses-s1Bx1" "ses-s1Bx1" "ses-s1Bx1" "ses-s1Bx1")

declare -a S1A_SUBS=("sub-cIIs00" "sub-cIIs01" "sub-cIIs02" "sub-cIIs03" "sub-cIIs04")
declare -a S1A_SES=("ses-s1Ax1" "ses-s1Ax1" "ses-s1Ax1" "ses-s1Ax1" "ses-s1Ax1")

declare -a S3_SUBS=("sub-cIIs00" "sub-cIIs01" "sub-cIIs02" "sub-cIIs03" "sub-cIIs04")
declare -a S3_SES=("ses-s3x1" "ses-s3x1" "ses-s3x1" "ses-s3x1" "ses-s3x1")

cd "$PROJECT_DIR"

# =============================================================================
# Step 1: Create covariate-adjusted template
# =============================================================================
echo "============================================================"
echo "Step 1: Creating covariate-adjusted RISH template"
echo "============================================================"

if [ -f "$TEMPLATE_COV_DIR/template_rish_l0.mif" ] && [ -f "$TEMPLATE_COV_DIR/covariate_model/covariate_model.json" ]; then
    echo "  [skip] Template already exists"
else
    mkdir -p "$TEMPLATE_COV_DIR"

    $SYS_PYTHON -m src.cli.main create-template \
        --reference-list "$SCRIPT_DIR/s1b_fods.txt" \
        --mask-list "$SCRIPT_DIR/s1b_masks.txt" \
        --output "$TEMPLATE_COV_DIR" \
        --lmax 8 \
        --nthreads "$NTHREADS" \
        --participants "$PARTICIPANTS" \
        --covariates age,sex \
        --subject-column participant_id

    echo ""
    echo "  Verifying outputs:"
    echo "    Template RISH files:"
    ls "$TEMPLATE_COV_DIR"/template_rish_l*.mif 2>/dev/null || echo "    ERROR: No template files"
    echo "    Covariate model:"
    ls "$TEMPLATE_COV_DIR"/covariate_model/covariate_model.json 2>/dev/null || echo "    ERROR: No model"
    echo "    Beta maps:"
    ls "$TEMPLATE_COV_DIR"/covariate_model/beta_*.mif 2>/dev/null || echo "    ERROR: No beta maps"
fi

echo ""

# =============================================================================
# Step 2: Harmonize target sites with covariate adjustment
# =============================================================================
echo "============================================================"
echo "Step 2: Harmonizing target sites with covariates"
echo "============================================================"

# Subject ages and sexes (for --subject-covariates)
declare -A SUBJ_AGE=( ["sub-cIIs00"]=37 ["sub-cIIs01"]=27 ["sub-cIIs02"]=45 ["sub-cIIs03"]=47 ["sub-cIIs04"]=36 )
declare -A SUBJ_SEX=( ["sub-cIIs00"]=M ["sub-cIIs01"]=M ["sub-cIIs02"]=F ["sub-cIIs03"]=F ["sub-cIIs04"]=M )

# Harmonize s1A subjects
for i in "${!S1A_SUBS[@]}"; do
    sub="${S1A_SUBS[$i]}"
    ses="${S1A_SES[$i]}"
    out="$HARM_COV_DIR/${sub}/${ses}"

    if [ -f "$out/sh_harmonized.mif" ]; then
        echo "  [skip] $sub $ses — harmonized"
        continue
    fi

    age="${SUBJ_AGE[$sub]}"
    sex="${SUBJ_SEX[$sub]}"

    echo "  Harmonizing $sub $ses (age=$age, sex=$sex) ..."
    mkdir -p "$out"

    $SYS_PYTHON -m src.cli.main harmonize \
        --target "$PROC_DIR/${sub}/${ses}/fod_reg.mif" \
        --mask "$PROC_DIR/group_mask.mif" \
        --template "$TEMPLATE_COV_DIR" \
        --output "$out" \
        --nthreads "$NTHREADS" \
        --subject-covariates "age=${age},sex=${sex}"
done

# Harmonize s3 subjects
for i in "${!S3_SUBS[@]}"; do
    sub="${S3_SUBS[$i]}"
    ses="${S3_SES[$i]}"
    out="$HARM_COV_DIR/${sub}/${ses}"

    if [ -f "$out/sh_harmonized.mif" ]; then
        echo "  [skip] $sub $ses — harmonized"
        continue
    fi

    age="${SUBJ_AGE[$sub]}"
    sex="${SUBJ_SEX[$sub]}"

    echo "  Harmonizing $sub $ses (age=$age, sex=$sex) ..."
    mkdir -p "$out"

    $SYS_PYTHON -m src.cli.main harmonize \
        --target "$PROC_DIR/${sub}/${ses}/fod_reg.mif" \
        --mask "$PROC_DIR/group_mask.mif" \
        --template "$TEMPLATE_COV_DIR" \
        --output "$out" \
        --nthreads "$NTHREADS" \
        --subject-covariates "age=${age},sex=${sex}"
done

echo ""

# =============================================================================
# Step 3: Extract post-harmonization RISH l0
# =============================================================================
echo "============================================================"
echo "Step 3: Extracting RISH l0 from covariate-harmonized FODs"
echo "============================================================"

# Reference site (s1B): use pre-harmonization l0 (unchanged)
for i in "${!S1B_SUBS[@]}"; do
    sub="${S1B_SUBS[$i]}"
    ses="${S1B_SES[$i]}"
    src_rish="$PROC_DIR/${sub}/${ses}/rish_l0_pre.mif"
    dst_rish="$PROC_DIR/${sub}/${ses}/rish_l0_postcov.mif"
    if [ ! -f "$dst_rish" ]; then
        cp "$src_rish" "$dst_rish"
    fi
done

# Target sites: extract l0 from covariate-harmonized FODs
for subs_var in S1A S3; do
    eval "subs=(\"\${${subs_var}_SUBS[@]}\")"
    eval "sess=(\"\${${subs_var}_SES[@]}\")"
    for i in "${!subs[@]}"; do
        sub="${subs[$i]}"
        ses="${sess[$i]}"
        harm_fod="$HARM_COV_DIR/${sub}/${ses}/sh_harmonized.mif"
        dst_rish="$PROC_DIR/${sub}/${ses}/rish_l0_postcov.mif"

        if [ -f "$dst_rish" ]; then
            echo "  [skip] $sub $ses"
            continue
        fi

        echo "  Extracting l0 for $sub $ses ..."
        mrconvert "$harm_fod" -coord 3 0 "$dst_rish" -force -quiet
    done
done

echo ""

# =============================================================================
# Step 4: Run site-effect test on covariate-harmonized data
# =============================================================================
echo "============================================================"
echo "Step 4: Site-effect test (post covariate-harmonization)"
echo "============================================================"

CSV="$SCRIPT_DIR/sites_postcov.csv"
echo "subject,site,image_path" > "$CSV"

declare -a ALL_SUBS=("${S1A_SUBS[@]}" "${S1B_SUBS[@]}" "${S3_SUBS[@]}")
declare -a ALL_SES=("${S1A_SES[@]}" "${S1B_SES[@]}" "${S3_SES[@]}")
declare -a ALL_SITES=()
for _ in "${S1A_SUBS[@]}"; do ALL_SITES+=("s1A"); done
for _ in "${S1B_SUBS[@]}"; do ALL_SITES+=("s1B"); done
for _ in "${S3_SUBS[@]}"; do ALL_SITES+=("s3"); done

for i in "${!ALL_SUBS[@]}"; do
    sub="${ALL_SUBS[$i]}"
    ses="${ALL_SES[$i]}"
    site="${ALL_SITES[$i]}"
    rish="$PROC_DIR/${sub}/${ses}/rish_l0_postcov.mif"
    echo "${sub}_${ses},${site},${rish}" >> "$CSV"
done

mkdir -p "$RESULTS_DIR/postcov"

$SYS_PYTHON -m src.cli.main site-effect \
    --site-list "$CSV" \
    --mask "$PROC_DIR/group_mask.mif" \
    --output "$RESULTS_DIR/postcov" \
    --n-permutations "$N_PERMS" \
    --seed 42

echo ""

# =============================================================================
# Step 5: Compare all three: pre / post-baseline / post-covariate
# =============================================================================
echo "============================================================"
echo "Step 5: Comparison — pre vs baseline-post vs covariate-post"
echo "============================================================"

$SYS_PYTHON -c "
import json
from pathlib import Path

pre = json.loads(Path('$RESULTS_DIR/pre/summary.json').read_text())
post = json.loads(Path('$RESULTS_DIR/post/summary.json').read_text())
postcov = json.loads(Path('$RESULTS_DIR/postcov/summary.json').read_text())

print('                         % Sig Voxels   Mean Effect Size')
print('                         -----------   ----------------')
print(f'  Pre-harmonization:      {pre[\"percent_significant_permutation\"]:>8.1f}%      {pre[\"mean_effect_size\"]:.4f}')
print(f'  Post (baseline):        {post[\"percent_significant_permutation\"]:>8.1f}%      {post[\"mean_effect_size\"]:.4f}')
print(f'  Post (covariate-adj):   {postcov[\"percent_significant_permutation\"]:>8.1f}%      {postcov[\"mean_effect_size\"]:.4f}')
print()

# Since same subjects on all scanners, covariate adjustment should give
# very similar results to baseline (no demographic confound to remove)
sig_diff = abs(post['percent_significant_permutation'] - postcov['percent_significant_permutation'])
es_diff = abs(post['mean_effect_size'] - postcov['mean_effect_size'])

print(f'  Baseline vs covariate-adjusted difference:')
print(f'    Sig voxels: {sig_diff:.1f} percentage points')
print(f'    Effect size: {es_diff:.4f}')
print()

if postcov['percent_significant_permutation'] < 5.0:
    print('  RESULT: Covariate-adjusted harmonization SUCCESSFUL')
    print('          (significant voxels < 5%)')
else:
    print('  RESULT: Covariate-adjusted harmonization needs review')
    print(f'          (significant voxels: {postcov[\"percent_significant_permutation\"]:.1f}%)')

print()
print('  NOTE: MASiVar is a within-subject design (same 5 subjects on all scanners),')
print('  so covariate adjustment should produce similar results to baseline.')
print('  The key verification is that the pipeline runs end-to-end and produces')
print('  valid covariate model files (beta maps, covariate_model.json).')

# Save comparison
comparison = {
    'pre_harmonization': {
        'percent_significant': pre['percent_significant_permutation'],
        'mean_effect_size': pre['mean_effect_size'],
    },
    'post_baseline': {
        'percent_significant': post['percent_significant_permutation'],
        'mean_effect_size': post['mean_effect_size'],
    },
    'post_covariate_adjusted': {
        'percent_significant': postcov['percent_significant_permutation'],
        'mean_effect_size': postcov['mean_effect_size'],
    },
    'baseline_vs_covariate_diff': {
        'sig_voxels_diff_pct': round(sig_diff, 2),
        'effect_size_diff': round(es_diff, 5),
    },
}
Path('$RESULTS_DIR/postcov/comparison.json').write_text(json.dumps(comparison, indent=2))
print()
print(f'  Full comparison saved to: $RESULTS_DIR/postcov/comparison.json')
"

echo ""
echo "============================================================"
echo "Covariate pipeline test complete"
echo "============================================================"

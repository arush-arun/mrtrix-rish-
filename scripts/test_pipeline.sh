#!/bin/bash
# Test pipeline on UKB sample data
# Tests: sub-002 and sub-003 from different sites

set -e

# Add MRtrix3 to PATH
export PATH="/home/uqahonne/mrtrix3/bin:$PATH"

# Create python symlink if needed (MRtrix3 scripts use 'python')
if ! command -v python &> /dev/null; then
    mkdir -p /tmp/pybin
    ln -sf $(which python3) /tmp/pybin/python
    export PATH="/tmp/pybin:$PATH"
fi

BIDS_DIR="/home/uqahonne/uq/ukb/data/sample_cud_dmri"
OUTPUT_DIR="/home/uqahonne/uq/nif/mrtrix-rish/test_output"
NTHREADS=4

mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "MRTRIX-RISH Pipeline Test"
echo "========================================"
echo ""

# Process sub-002 ses-1
process_subject() {
    local SUBJ=$1
    local SES=$2
    local SITE=$3
    
    echo "Processing $SUBJ / $SES (Site: $SITE)..."
    
    local DWI_DIR="$BIDS_DIR/$SUBJ/$SES/dwi"
    local DWI_NII="$DWI_DIR/${SUBJ}_${SES}_dir-AP_dwi.nii.gz"
    local BVAL="$DWI_DIR/${SUBJ}_${SES}_dir-AP_dwi.bval"
    local BVEC="$DWI_DIR/${SUBJ}_${SES}_dir-AP_dwi.bvec"
    
    local OUT="$OUTPUT_DIR/$SUBJ/$SES"
    mkdir -p "$OUT"
    
    # Step 1: Convert to MIF
    echo "  [1/5] Converting to MIF..."
    local MIF="$OUT/dwi.mif"
    if [ ! -f "$MIF" ]; then
        mrconvert "$DWI_NII" -fslgrad "$BVEC" "$BVAL" "$MIF" -force -nthreads $NTHREADS
    fi
    
    # Step 2: Create mask
    echo "  [2/5] Creating brain mask..."
    local MASK="$OUT/mask.mif"
    if [ ! -f "$MASK" ]; then
        dwi2mask "$MIF" "$MASK" -force -nthreads $NTHREADS
    fi
    
    # Step 3: Check shells
    echo "  [3/5] Detecting shells..."
    mrinfo -shell_bvalues "$MIF"
    mrinfo -shell_sizes "$MIF"
    
    # Step 4: Estimate response (single-shell -> tournier)
    echo "  [4/5] Estimating response function..."
    local RESPONSE="$OUT/response_wm.txt"
    if [ ! -f "$RESPONSE" ]; then
        dwi2response tournier "$MIF" "$RESPONSE" -mask "$MASK" -force -nthreads $NTHREADS
    fi
    
    # Step 5: Compute FOD (single-shell CSD)
    echo "  [5/5] Computing FOD (CSD)..."
    local FOD="$OUT/wm_fod.mif"
    if [ ! -f "$FOD" ]; then
        dwi2fod csd "$MIF" "$RESPONSE" "$FOD" -mask "$MASK" -lmax 8 -force -nthreads $NTHREADS
    fi
    
    echo "  ✓ Done: $FOD"
    echo ""
}

# Process both subjects
process_subject "sub-002" "ses-1" "Site_A"
process_subject "sub-003" "ses-1" "Site_B"

echo "========================================"
echo "FODs computed. Now extracting RISH..."
echo "========================================"

# Extract RISH features
extract_rish() {
    local SUBJ=$1
    local SES=$2
    
    local FOD="$OUTPUT_DIR/$SUBJ/$SES/wm_fod.mif"
    local MASK="$OUTPUT_DIR/$SUBJ/$SES/mask.mif"
    local RISH_DIR="$OUTPUT_DIR/$SUBJ/$SES/rish"
    
    mkdir -p "$RISH_DIR"
    
    echo "Extracting RISH for $SUBJ..."
    
    # Get lmax from FOD
    NVOLS=$(mrinfo -size "$FOD" | awk '{print $4}')
    echo "  FOD has $NVOLS volumes"
    
    # For lmax=8: volumes 0-44 (45 total)
    # l=0: vol 0 (1 coeff)
    # l=2: vol 1-5 (5 coeffs)
    # l=4: vol 6-14 (9 coeffs)
    # l=6: vol 15-27 (13 coeffs)
    # l=8: vol 28-44 (17 coeffs)
    
    # θ_0 = c_00 (just volume 0)
    echo "  Extracting θ_0..."
    mrconvert "$FOD" -coord 3 0 "$RISH_DIR/rish_l0.mif" -force -quiet
    mrcalc "$RISH_DIR/rish_l0.mif" "$MASK" -mult "$RISH_DIR/rish_l0.mif" -force -quiet
    
    # θ_2 = sqrt(sum(c_2m^2)) for m=-2..2 (volumes 1-5)
    echo "  Extracting θ_2..."
    mrconvert "$FOD" -coord 3 1:5 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow - -quiet | mrcalc - "$MASK" -mult "$RISH_DIR/rish_l2.mif" -force -quiet
    
    # θ_4 (volumes 6-14)
    echo "  Extracting θ_4..."
    mrconvert "$FOD" -coord 3 6:14 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow - -quiet | mrcalc - "$MASK" -mult "$RISH_DIR/rish_l4.mif" -force -quiet
    
    # θ_6 (volumes 15-27)
    echo "  Extracting θ_6..."
    mrconvert "$FOD" -coord 3 15:27 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow - -quiet | mrcalc - "$MASK" -mult "$RISH_DIR/rish_l6.mif" -force -quiet
    
    # θ_8 (volumes 28-44)
    echo "  Extracting θ_8..."
    mrconvert "$FOD" -coord 3 28:44 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow - -quiet | mrcalc - "$MASK" -mult "$RISH_DIR/rish_l8.mif" -force -quiet
    
    echo "  ✓ RISH features: $RISH_DIR/"
    ls "$RISH_DIR/"
    echo ""
}

extract_rish "sub-002" "ses-1"
extract_rish "sub-003" "ses-1"

echo "========================================"
echo "Computing RISH statistics..."
echo "========================================"

for SUBJ in sub-002 sub-003; do
    echo "$SUBJ RISH stats:"
    for L in 0 2 4 6 8; do
        MEAN=$(mrstats "$OUTPUT_DIR/$SUBJ/ses-1/rish/rish_l${L}.mif" -mask "$OUTPUT_DIR/$SUBJ/ses-1/mask.mif" -output mean 2>/dev/null)
        echo "  θ_$L mean: $MEAN"
    done
    echo ""
done

echo "========================================"
echo "✓ Pipeline test complete!"
echo "========================================"
echo ""
echo "Output: $OUTPUT_DIR"
echo ""
echo "Next steps:"
echo "  1. Compare RISH features between sites"
echo "  2. Compute scale maps for harmonization"

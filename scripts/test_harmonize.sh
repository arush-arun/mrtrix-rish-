#!/bin/bash
# Test RISH Harmonization
# Harmonize sub-003 (Site B) to sub-002 (Site A) as reference

set -e

export PATH="/home/uqahonne/mrtrix3/bin:$PATH"
if ! command -v python &> /dev/null; then
    mkdir -p /tmp/pybin
    ln -sf $(which python3) /tmp/pybin/python
    export PATH="/tmp/pybin:$PATH"
fi

OUTPUT_DIR="/home/uqahonne/uq/nif/mrtrix-rish/test_output"
HARMONIZE_DIR="$OUTPUT_DIR/harmonization"
NTHREADS=4

mkdir -p "$HARMONIZE_DIR"

echo "========================================"
echo "RISH HARMONIZATION"
echo "========================================"
echo ""
echo "Reference: sub-002 (Site A)"
echo "Target:    sub-003 (Site B)"
echo ""

# Paths
REF_RISH="$OUTPUT_DIR/sub-002/ses-1/rish"
TAR_RISH="$OUTPUT_DIR/sub-003/ses-1/rish"
TAR_FOD="$OUTPUT_DIR/sub-003/ses-1/wm_fod.mif"
TAR_MASK="$OUTPUT_DIR/sub-003/ses-1/mask.mif"

# Step 1: Compute scale maps for each SH order
echo "[1/3] Computing scale maps..."
SCALE_DIR="$HARMONIZE_DIR/scale_maps"
mkdir -p "$SCALE_DIR"

for L in 0 2 4 6 8; do
    echo "  Computing scale_l$L..."
    REF="$REF_RISH/rish_l${L}.mif"
    TAR="$TAR_RISH/rish_l${L}.mif"
    SCALE="$SCALE_DIR/scale_l${L}.mif"
    
    # scale = reference / target (with smoothing and clipping)
    # Step 1: Threshold target to avoid division by zero
    mrcalc "$TAR" 0.01 -max /tmp/tar_thresh.mif -force -quiet
    
    # Step 2: Compute ratio
    mrcalc "$REF" /tmp/tar_thresh.mif -div /tmp/ratio.mif -force -quiet
    
    # Step 3: Apply mask
    mrcalc /tmp/ratio.mif "$TAR_MASK" -mult /tmp/ratio_masked.mif -force -quiet
    
    # Step 4: Smooth (sigma = 3mm / 2.355 ≈ 1.27)
    mrfilter /tmp/ratio_masked.mif smooth -stdev 1.27 /tmp/ratio_smooth.mif -force -quiet
    
    # Step 5: Clip to [0.5, 2.0]
    mrcalc /tmp/ratio_smooth.mif 0.5 -max 2.0 -min /tmp/ratio_clipped.mif -force -quiet
    
    # Step 6: Set non-brain to 1.0
    mrcalc "$TAR_MASK" /tmp/ratio_clipped.mif -mult \
           "$TAR_MASK" 1 -sub -neg -add \
           "$SCALE" -force -quiet
    
    # Stats
    MEAN=$(mrstats "$SCALE" -mask "$TAR_MASK" -output mean 2>/dev/null)
    echo "    scale_l$L mean: $MEAN"
done

echo ""

# Step 2: Apply harmonization to FOD
echo "[2/3] Applying harmonization to FOD..."
HARM_FOD="$HARMONIZE_DIR/wm_fod_harmonized.mif"

# Extract and scale each SH order, then concatenate
# l=0: vol 0 (1 coeff)
# l=2: vol 1-5 (5 coeffs)
# l=4: vol 6-14 (9 coeffs)
# l=6: vol 15-27 (13 coeffs)
# l=8: vol 28-44 (17 coeffs)

echo "  Scaling l=0..."
mrconvert "$TAR_FOD" -coord 3 0 - -quiet | mrcalc - "$SCALE_DIR/scale_l0.mif" -mult /tmp/sh_l0.mif -force -quiet

echo "  Scaling l=2..."
mrconvert "$TAR_FOD" -coord 3 1:5 - -quiet | mrcalc - "$SCALE_DIR/scale_l2.mif" -mult /tmp/sh_l2.mif -force -quiet

echo "  Scaling l=4..."
mrconvert "$TAR_FOD" -coord 3 6:14 - -quiet | mrcalc - "$SCALE_DIR/scale_l4.mif" -mult /tmp/sh_l4.mif -force -quiet

echo "  Scaling l=6..."
mrconvert "$TAR_FOD" -coord 3 15:27 - -quiet | mrcalc - "$SCALE_DIR/scale_l6.mif" -mult /tmp/sh_l6.mif -force -quiet

echo "  Scaling l=8..."
mrconvert "$TAR_FOD" -coord 3 28:44 - -quiet | mrcalc - "$SCALE_DIR/scale_l8.mif" -mult /tmp/sh_l8.mif -force -quiet

echo "  Concatenating..."
mrcat /tmp/sh_l0.mif /tmp/sh_l2.mif /tmp/sh_l4.mif /tmp/sh_l6.mif /tmp/sh_l8.mif \
      -axis 3 "$HARM_FOD" -force -quiet

echo "  ✓ Harmonized FOD: $HARM_FOD"
echo ""

# Step 3: Extract RISH from harmonized FOD
echo "[3/3] Extracting RISH from harmonized FOD..."
HARM_RISH="$HARMONIZE_DIR/rish_harmonized"
mkdir -p "$HARM_RISH"

mrconvert "$HARM_FOD" -coord 3 0 "$HARM_RISH/rish_l0.mif" -force -quiet
mrconvert "$HARM_FOD" -coord 3 1:5 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow "$HARM_RISH/rish_l2.mif" -force -quiet
mrconvert "$HARM_FOD" -coord 3 6:14 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow "$HARM_RISH/rish_l4.mif" -force -quiet
mrconvert "$HARM_FOD" -coord 3 15:27 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow "$HARM_RISH/rish_l6.mif" -force -quiet
mrconvert "$HARM_FOD" -coord 3 28:44 - -quiet | mrcalc - 2 -pow - -quiet | mrmath - sum -axis 3 - -quiet | mrcalc - 0.5 -pow "$HARM_RISH/rish_l8.mif" -force -quiet

echo ""
echo "========================================"
echo "COMPARISON: Before vs After Harmonization"
echo "========================================"
echo ""
echo "                    Reference     Original      Harmonized    Δ Original   Δ Harmonized"
echo "                    (sub-002)     (sub-003)     (sub-003)     vs Ref       vs Ref"
echo "--------------------------------------------------------------------------------"

for L in 0 2 4 6 8; do
    REF_MEAN=$(mrstats "$REF_RISH/rish_l${L}.mif" -mask "$OUTPUT_DIR/sub-002/ses-1/mask.mif" -output mean 2>/dev/null)
    TAR_MEAN=$(mrstats "$TAR_RISH/rish_l${L}.mif" -mask "$TAR_MASK" -output mean 2>/dev/null)
    HARM_MEAN=$(mrstats "$HARM_RISH/rish_l${L}.mif" -mask "$TAR_MASK" -output mean 2>/dev/null)
    
    # Calculate percentage differences
    DELTA_ORIG=$(echo "scale=1; ($TAR_MEAN - $REF_MEAN) / $REF_MEAN * 100" | bc 2>/dev/null || echo "N/A")
    DELTA_HARM=$(echo "scale=1; ($HARM_MEAN - $REF_MEAN) / $REF_MEAN * 100" | bc 2>/dev/null || echo "N/A")
    
    printf "θ_%d                 %-12.4f  %-12.4f  %-12.4f  %+6.1f%%       %+6.1f%%\n" \
           $L $REF_MEAN $TAR_MEAN $HARM_MEAN $DELTA_ORIG $DELTA_HARM
done

echo ""
echo "========================================"
echo "✓ Harmonization complete!"
echo "========================================"
echo ""
echo "Output files:"
echo "  Scale maps:     $SCALE_DIR/"
echo "  Harmonized FOD: $HARM_FOD"
echo "  Harmonized RISH: $HARM_RISH/"
echo ""
echo "The harmonized sub-003 RISH features should now be closer to sub-002."

# Cleanup
rm -f /tmp/tar_thresh.mif /tmp/ratio.mif /tmp/ratio_masked.mif /tmp/ratio_smooth.mif /tmp/ratio_clipped.mif
rm -f /tmp/sh_l0.mif /tmp/sh_l2.mif /tmp/sh_l4.mif /tmp/sh_l6.mif /tmp/sh_l8.mif

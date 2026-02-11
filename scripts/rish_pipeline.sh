#!/bin/bash

# RISH Harmonization Pipeline
# Author: Arush Arun / Axon
# Complete preprocessing + RISH harmonization for multi-site dMRI
# Based on MRtrix3_FBA workflow structure

set -e

# Add MRtrix3 to PATH
export PATH="/home/uqahonne/mrtrix3/bin:$PATH"

# Python symlink for MRtrix3 scripts
if ! command -v python &> /dev/null; then
    mkdir -p /tmp/pybin
    ln -sf $(which python3) /tmp/pybin/python
    export PATH="/tmp/pybin:$PATH"
fi

# Default config file
CONFIG_FILE="${CONFIG_FILE:-./scripts/rish_config.conf}"

#=============================================================================
# CONFIGURATION
#=============================================================================

load_config() {
    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "ERROR: Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
    
    echo "Loading configuration from: $CONFIG_FILE"
    
    while IFS= read -r line; do
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
        export "$line" 2>/dev/null || true
    done < <(grep -E '^[A-Z_].*=' "$CONFIG_FILE")
    
    # Set derived paths
    if [[ "$DERIVATIVES_PATH" =~ ^/ ]]; then
        DER_PATH="$DERIVATIVES_PATH"
    else
        DER_PATH="${INPUT_PATH}/${DERIVATIVES_PATH}"
    fi
    
    echo "✓ Configuration loaded"
}

show_config() {
    echo "=============================================="
    echo "RISH HARMONIZATION PIPELINE CONFIGURATION"
    echo "=============================================="
    echo "Input path:        $INPUT_PATH"
    echo "Derivatives:       $DER_PATH"
    echo "Session:           $SESSION"
    echo "PE directions:     ${MAIN_PE_DIR} / ${BLIPPED_PE_DIR}"
    echo "Voxel size:        $VOXEL_SIZE"
    echo "lmax:              $LMAX"
    echo "Response algo:     $RESPONSE_ALGORITHM"
    echo "Reference:         $REFERENCE_SUBJECTS"
    echo "Targets:           $TARGET_SUBJECTS"
    echo "Smoothing FWHM:    $SMOOTHING_FWHM"
    echo "Scale range:       [$SCALE_MIN, $SCALE_MAX]"
    echo "Threads:           $NTHREADS"
    echo "=============================================="
}

#=============================================================================
# CHECKPOINT FUNCTIONS
#=============================================================================

create_checkpoint() {
    [[ "$ENABLE_CHECKPOINTS" != "true" ]] && return 0
    local subject="$1"
    local step="$2"
    mkdir -p "$CHECKPOINT_DIR"
    touch "${CHECKPOINT_DIR}/${subject}_${step}.done"
    echo "  ✓ Checkpoint: ${subject}/${step}"
}

check_checkpoint() {
    [[ "$ENABLE_CHECKPOINTS" != "true" ]] && return 1
    local subject="$1"
    local step="$2"
    [[ -f "${CHECKPOINT_DIR}/${subject}_${step}.done" ]]
}

#=============================================================================
# FILE PATH HELPERS
#=============================================================================

get_bids_files() {
    local subject="$1"
    local pe_dir="$2"
    local subject_dir="${INPUT_PATH}/${subject}/${SESSION}/dwi"
    
    echo "${subject_dir}/${subject}_${SESSION}_dir-${pe_dir}_dwi.nii.gz"
    echo "${subject_dir}/${subject}_${SESSION}_dir-${pe_dir}_dwi.bval"
    echo "${subject_dir}/${subject}_${SESSION}_dir-${pe_dir}_dwi.bvec"
}

#=============================================================================
# PREPROCESSING
#=============================================================================

preprocess_subject() {
    local subject="$1"
    local out_dir="${DER_PATH}/${subject}/${SESSION}/dwi"
    
    echo ""
    echo "############################################################"
    echo "PREPROCESSING: $subject"
    echo "############################################################"
    
    if check_checkpoint "$subject" "preprocess"; then
        echo "⏭️  Already completed, skipping..."
        return 0
    fi
    
    mkdir -p "$out_dir"
    
    # Get file paths
    local main_files=($(get_bids_files "$subject" "$MAIN_PE_DIR"))
    local blip_files=($(get_bids_files "$subject" "$BLIPPED_PE_DIR"))
    
    local main_dwi="${main_files[0]}"
    local main_bval="${main_files[1]}"
    local main_bvec="${main_files[2]}"
    local blip_dwi="${blip_files[0]}"
    local blip_bval="${blip_files[1]}"
    local blip_bvec="${blip_files[2]}"
    
    echo "  Main DWI: $main_dwi"
    echo "  Blip DWI: $blip_dwi"
    
    # Step 1: Convert to MIF
    echo "  [1/7] Converting to MIF..."
    mrconvert "$main_dwi" "$out_dir/dwi_raw.mif" \
              -fslgrad "$main_bvec" "$main_bval" \
              -nthreads $NTHREADS -force -quiet
    
    # Step 2: Denoise (need an option of adding denosing from DIPY)
    echo "  [2/7] Denoising..."
    dwidenoise "$out_dir/dwi_raw.mif" "$out_dir/dwi_denoised.mif" \
               -noise "$out_dir/noise.mif" \
               -nthreads $NTHREADS -force -quiet
    
    # Calculate residual and SNR
    mrcalc "$out_dir/dwi_raw.mif" "$out_dir/dwi_denoised.mif" \
           -subtract "$out_dir/residual_denoise.mif" -force -quiet
    
    # Step 3: Gibbs unringing
    echo "  [3/7] Gibbs unringing..."
    mrdegibbs "$out_dir/dwi_denoised.mif" "$out_dir/dwi_unring.mif" \
              -axes 0,1 -nthreads $NTHREADS -force -quiet
    
    # Step 4: Prepare for distortion correction
    echo "  [4/7] Preparing distortion correction..."
    
    # Convert blipped data
    mrconvert "$blip_dwi" "$out_dir/dwi_blip.mif" \
              -fslgrad "$blip_bvec" "$blip_bval" \
              -nthreads $NTHREADS -force -quiet
    
    # Extract and average b0s
    dwiextract "$out_dir/dwi_unring.mif" -bzero - -quiet | \
        mrmath - mean "$out_dir/b0_main.mif" -axis 3 -force -quiet
    
    dwiextract "$out_dir/dwi_blip.mif" -bzero - -quiet | \
        mrmath - mean "$out_dir/b0_blip.mif" -axis 3 -force -quiet
    
    mrcat "$out_dir/b0_main.mif" "$out_dir/b0_blip.mif" \
          "$out_dir/b0_pair.mif" -axis 3 -force -quiet
    
    # Step 5: Motion + Distortion correction
    echo "  [5/7] Running dwifslpreproc (this takes a while)..."
    dwifslpreproc "$out_dir/dwi_unring.mif" "$out_dir/dwi_preproc.mif" \
                  -pe_dir "$PHASE_ENCODING_DIR" -rpe_pair \
                  -se_epi "$out_dir/b0_pair.mif" \
                  -eddy_options "$EDDY_OPTIONS" \
                  -eddyqc_all "$out_dir/eddy_qc" \
                  -nthreads $NTHREADS -force -quiet
    
    # Step 6: Bias field correction
    echo "  [6/7] Bias field correction..."
    dwibiascorrect ants "$out_dir/dwi_preproc.mif" "$out_dir/dwi_unbiased.mif" \
                   -bias "$out_dir/bias.mif" \
                   -nthreads $NTHREADS -force -quiet
    
    # Step 7: Create brain mask
    echo "  [7/7] Creating brain mask..."
    dwi2mask "$out_dir/dwi_unbiased.mif" "$out_dir/mask.mif" \
             -nthreads $NTHREADS -force -quiet
    
    # Optional: Upsample
    if [[ "$VOXEL_SIZE" != "native" ]]; then
        echo "  Upsampling to ${VOXEL_SIZE}mm..."
        mrgrid "$out_dir/dwi_unbiased.mif" regrid \
               -vox "$VOXEL_SIZE" "$out_dir/dwi_final.mif" \
               -nthreads $NTHREADS -force -quiet
        mrgrid "$out_dir/mask.mif" regrid \
               -vox "$VOXEL_SIZE" "$out_dir/mask_upsampled.mif" \
               -interp nearest -nthreads $NTHREADS -force -quiet
    else
        ln -sf dwi_unbiased.mif "$out_dir/dwi_final.mif"
        ln -sf mask.mif "$out_dir/mask_upsampled.mif"
    fi
    
    # Calculate SNR
    local SNR=$(dwiextract "$out_dir/dwi_final.mif" -no_bzero -singleshell - -quiet 2>/dev/null | \
                mrcalc - "$out_dir/noise.mif" -div - -quiet 2>/dev/null | \
                mrstats - -mask "$out_dir/mask.mif" -output mean -allvolumes 2>/dev/null || echo "N/A")
    echo "  SNR: $SNR"
    
    create_checkpoint "$subject" "preprocess"
    echo "  ✓ Preprocessing complete: $subject"
}

#=============================================================================
# FOD COMPUTATION
#=============================================================================

compute_fod_subject() {
    local subject="$1"
    local out_dir="${DER_PATH}/${subject}/${SESSION}/dwi"
    
    echo ""
    echo "############################################################"
    echo "FOD COMPUTATION: $subject"
    echo "############################################################"
    
    if check_checkpoint "$subject" "fod"; then
        echo "⏭️  Already completed, skipping..."
        return 0
    fi
    
    local dwi="$out_dir/dwi_final.mif"
    local mask="$out_dir/mask.mif"
    
    # Detect shell structure
    echo "  Detecting shell structure..."
    local shells=$(mrinfo -shell_bvalues "$dwi" 2>/dev/null)
    local n_shells=$(echo "$shells" | wc -w)
    
    echo "  Shells: $shells (n=$n_shells)"
    
    # Determine response algorithm
    local resp_algo="$RESPONSE_ALGORITHM"
    if [[ "$resp_algo" == "auto" ]]; then
        if [[ $n_shells -gt 2 ]]; then
            resp_algo="dhollander"
        else
            resp_algo="tournier"
        fi
    fi
    echo "  Response algorithm: $resp_algo"
    
    # Estimate response function
    echo "  [1/2] Estimating response function..."
    if [[ "$resp_algo" == "dhollander" ]]; then
        dwi2response dhollander "$dwi" \
                     "$out_dir/response_wm.txt" \
                     "$out_dir/response_gm.txt" \
                     "$out_dir/response_csf.txt" \
                     -voxels "$out_dir/response_voxels.mif" \
                     -mask "$mask" -nthreads $NTHREADS -force -quiet
        
        echo "  [2/2] Computing FOD (MSMT-CSD)..."
        dwi2fod msmt_csd "$dwi" \
                "$out_dir/response_wm.txt" "$out_dir/wm_fod.mif" \
                "$out_dir/response_gm.txt" "$out_dir/gm.mif" \
                "$out_dir/response_csf.txt" "$out_dir/csf.mif" \
                -mask "$mask" -lmax $LMAX,$LMAX,$LMAX \
                -nthreads $NTHREADS -force -quiet
    else
        dwi2response tournier "$dwi" \
                     "$out_dir/response_wm.txt" \
                     -voxels "$out_dir/response_voxels.mif" \
                     -mask "$mask" -nthreads $NTHREADS -force -quiet
        
        echo "  [2/2] Computing FOD (CSD)..."
        dwi2fod csd "$dwi" \
                "$out_dir/response_wm.txt" "$out_dir/wm_fod.mif" \
                -mask "$mask" -lmax $LMAX \
                -nthreads $NTHREADS -force -quiet
    fi
    
    create_checkpoint "$subject" "fod"
    echo "  ✓ FOD complete: $subject"
}

#=============================================================================
# RISH EXTRACTION
#=============================================================================

extract_rish() {
    local subject="$1"
    local out_dir="${DER_PATH}/${subject}/${SESSION}/dwi"
    local rish_dir="$out_dir/rish"
    
    echo ""
    echo "############################################################"
    echo "RISH EXTRACTION: $subject"
    echo "############################################################"
    
    if check_checkpoint "$subject" "rish_extract"; then
        echo "⏭️  Already completed, skipping..."
        return 0
    fi
    
    mkdir -p "$rish_dir"
    
    local fod="$out_dir/wm_fod.mif"
    local mask="$out_dir/mask.mif"
    
    # Get number of volumes to determine lmax
    local nvols=$(mrinfo -size "$fod" | awk '{print $4}')
    echo "  FOD volumes: $nvols"
    
    # Extract RISH for each SH order
    # l=0: vol 0
    echo "  Extracting θ₀..."
    mrconvert "$fod" -coord 3 0 "$rish_dir/rish_l0.mif" -force -quiet
    mrcalc "$rish_dir/rish_l0.mif" "$mask" -mult "$rish_dir/rish_l0.mif" -force -quiet
    
    # l=2: vol 1:5
    echo "  Extracting θ₂..."
    mrconvert "$fod" -coord 3 1:5 - -quiet | \
        mrcalc - 2 -pow - -quiet | \
        mrmath - sum -axis 3 - -quiet | \
        mrcalc - 0.5 -pow - -quiet | \
        mrcalc - "$mask" -mult "$rish_dir/rish_l2.mif" -force -quiet
    
    # l=4: vol 6:14
    echo "  Extracting θ₄..."
    mrconvert "$fod" -coord 3 6:14 - -quiet | \
        mrcalc - 2 -pow - -quiet | \
        mrmath - sum -axis 3 - -quiet | \
        mrcalc - 0.5 -pow - -quiet | \
        mrcalc - "$mask" -mult "$rish_dir/rish_l4.mif" -force -quiet
    
    # l=6: vol 15:27
    echo "  Extracting θ₆..."
    mrconvert "$fod" -coord 3 15:27 - -quiet | \
        mrcalc - 2 -pow - -quiet | \
        mrmath - sum -axis 3 - -quiet | \
        mrcalc - 0.5 -pow - -quiet | \
        mrcalc - "$mask" -mult "$rish_dir/rish_l6.mif" -force -quiet
    
    # l=8: vol 28:44 (if lmax=8)
    if [[ $nvols -ge 45 ]]; then
        echo "  Extracting θ₈..."
        mrconvert "$fod" -coord 3 28:44 - -quiet | \
            mrcalc - 2 -pow - -quiet | \
            mrmath - sum -axis 3 - -quiet | \
            mrcalc - 0.5 -pow - -quiet | \
            mrcalc - "$mask" -mult "$rish_dir/rish_l8.mif" -force -quiet
    fi
    
    create_checkpoint "$subject" "rish_extract"
    echo "  ✓ RISH extraction complete: $subject"
}

#=============================================================================
# RISH HARMONIZATION
#=============================================================================

create_rish_template() {
    echo ""
    echo "############################################################"
    echo "CREATING RISH TEMPLATE FROM REFERENCE SUBJECTS"
    echo "############################################################"
    
    local template_dir="${DER_PATH}/template"
    mkdir -p "$template_dir"
    
    # Average RISH features across reference subjects
    for l in 0 2 4 6 8; do
        local rish_files=""
        for ref_subj in $REFERENCE_SUBJECTS; do
            local rish_file="${DER_PATH}/${ref_subj}/${SESSION}/dwi/rish/rish_l${l}.mif"
            if [[ -f "$rish_file" ]]; then
                rish_files="$rish_files $rish_file"
            fi
        done
        
        if [[ -n "$rish_files" ]]; then
            echo "  Averaging θ_${l} from reference subjects..."
            mrmath $rish_files mean "$template_dir/template_rish_l${l}.mif" \
                   -force -quiet -nthreads $NTHREADS
        fi
    done
    
    echo "  ✓ Template created: $template_dir"
}

harmonize_subject() {
    local subject="$1"
    local out_dir="${DER_PATH}/${subject}/${SESSION}/dwi"
    local harm_dir="$out_dir/harmonized"
    local template_dir="${DER_PATH}/template"
    
    echo ""
    echo "############################################################"
    echo "HARMONIZING: $subject"
    echo "############################################################"
    
    if check_checkpoint "$subject" "harmonize"; then
        echo "⏭️  Already completed, skipping..."
        return 0
    fi
    
    mkdir -p "$harm_dir/scale_maps"
    mkdir -p "$harm_dir/rish"
    
    local fod="$out_dir/wm_fod.mif"
    local mask="$out_dir/mask.mif"
    local rish_dir="$out_dir/rish"
    
    # Compute scale maps
    echo "  Computing scale maps..."
    local sigma=$(echo "scale=4; $SMOOTHING_FWHM / 2.355" | bc)
    
    for l in 0 2 4 6 8; do
        local ref="$template_dir/template_rish_l${l}.mif"
        local tar="$rish_dir/rish_l${l}.mif"
        local scale="$harm_dir/scale_maps/scale_l${l}.mif"
        
        if [[ -f "$ref" && -f "$tar" ]]; then
            echo "    scale_l${l}..."
            
            # Threshold target to avoid division by zero
            mrcalc "$tar" 0.01 -max /tmp/tar_thresh_$$.mif -force -quiet
            
            # Compute ratio
            mrcalc "$ref" /tmp/tar_thresh_$$.mif -div /tmp/ratio_$$.mif -force -quiet
            
            # Apply mask
            mrcalc /tmp/ratio_$$.mif "$mask" -mult /tmp/ratio_masked_$$.mif -force -quiet
            
            # Smooth
            mrfilter /tmp/ratio_masked_$$.mif smooth -stdev $sigma \
                     /tmp/ratio_smooth_$$.mif -force -quiet
            
            # Clip
            mrcalc /tmp/ratio_smooth_$$.mif $SCALE_MIN -max $SCALE_MAX -min \
                   /tmp/ratio_clipped_$$.mif -force -quiet
            
            # Set non-brain to 1.0
            mrcalc "$mask" /tmp/ratio_clipped_$$.mif -mult \
                   "$mask" 1 -sub -neg -add \
                   "$scale" -force -quiet
            
            rm -f /tmp/*_$$.mif
        fi
    done
    
    # Apply harmonization to FOD
    echo "  Applying harmonization to FOD..."
    
    mrconvert "$fod" -coord 3 0 - -quiet | \
        mrcalc - "$harm_dir/scale_maps/scale_l0.mif" -mult /tmp/sh_l0_$$.mif -force -quiet
    
    mrconvert "$fod" -coord 3 1:5 - -quiet | \
        mrcalc - "$harm_dir/scale_maps/scale_l2.mif" -mult /tmp/sh_l2_$$.mif -force -quiet
    
    mrconvert "$fod" -coord 3 6:14 - -quiet | \
        mrcalc - "$harm_dir/scale_maps/scale_l4.mif" -mult /tmp/sh_l4_$$.mif -force -quiet
    
    mrconvert "$fod" -coord 3 15:27 - -quiet | \
        mrcalc - "$harm_dir/scale_maps/scale_l6.mif" -mult /tmp/sh_l6_$$.mif -force -quiet
    
    mrconvert "$fod" -coord 3 28:44 - -quiet | \
        mrcalc - "$harm_dir/scale_maps/scale_l8.mif" -mult /tmp/sh_l8_$$.mif -force -quiet
    
    mrcat /tmp/sh_l0_$$.mif /tmp/sh_l2_$$.mif /tmp/sh_l4_$$.mif \
          /tmp/sh_l6_$$.mif /tmp/sh_l8_$$.mif \
          -axis 3 "$harm_dir/wm_fod_harmonized.mif" -force -quiet
    
    rm -f /tmp/sh_l*_$$.mif
    
    # Extract RISH from harmonized FOD
    echo "  Extracting harmonized RISH..."
    mrconvert "$harm_dir/wm_fod_harmonized.mif" -coord 3 0 \
              "$harm_dir/rish/rish_l0.mif" -force -quiet
    
    for l in 2 4 6 8; do
        case $l in
            2) coords="1:5" ;;
            4) coords="6:14" ;;
            6) coords="15:27" ;;
            8) coords="28:44" ;;
        esac
        
        mrconvert "$harm_dir/wm_fod_harmonized.mif" -coord 3 $coords - -quiet | \
            mrcalc - 2 -pow - -quiet | \
            mrmath - sum -axis 3 - -quiet | \
            mrcalc - 0.5 -pow "$harm_dir/rish/rish_l${l}.mif" -force -quiet
    done
    
    create_checkpoint "$subject" "harmonize"
    echo "  ✓ Harmonization complete: $subject"
}

#=============================================================================
# QC REPORT
#=============================================================================

generate_qc_report() {
    echo ""
    echo "############################################################"
    echo "GENERATING QC REPORT"
    echo "############################################################"
    
    local qc_dir="${DER_PATH}/qc"
    mkdir -p "$qc_dir"
    
    # Collect stats for all subjects
    local report_file="$qc_dir/harmonization_summary.csv"
    echo "subject,type,theta_0,theta_2,theta_4,theta_6,theta_8" > "$report_file"
    
    # Reference subjects
    for subj in $REFERENCE_SUBJECTS; do
        local rish_dir="${DER_PATH}/${subj}/${SESSION}/dwi/rish"
        local mask="${DER_PATH}/${subj}/${SESSION}/dwi/mask.mif"
        
        if [[ -d "$rish_dir" ]]; then
            local r0=$(mrstats "$rish_dir/rish_l0.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r2=$(mrstats "$rish_dir/rish_l2.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r4=$(mrstats "$rish_dir/rish_l4.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r6=$(mrstats "$rish_dir/rish_l6.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r8=$(mrstats "$rish_dir/rish_l8.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            echo "$subj,reference,$r0,$r2,$r4,$r6,$r8" >> "$report_file"
        fi
    done
    
    # Target subjects (original and harmonized)
    for subj in $TARGET_SUBJECTS; do
        local mask="${DER_PATH}/${subj}/${SESSION}/dwi/mask.mif"
        
        # Original
        local rish_dir="${DER_PATH}/${subj}/${SESSION}/dwi/rish"
        if [[ -d "$rish_dir" ]]; then
            local r0=$(mrstats "$rish_dir/rish_l0.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r2=$(mrstats "$rish_dir/rish_l2.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r4=$(mrstats "$rish_dir/rish_l4.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r6=$(mrstats "$rish_dir/rish_l6.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r8=$(mrstats "$rish_dir/rish_l8.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            echo "$subj,original,$r0,$r2,$r4,$r6,$r8" >> "$report_file"
        fi
        
        # Harmonized
        local harm_dir="${DER_PATH}/${subj}/${SESSION}/dwi/harmonized/rish"
        if [[ -d "$harm_dir" ]]; then
            local r0=$(mrstats "$harm_dir/rish_l0.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r2=$(mrstats "$harm_dir/rish_l2.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r4=$(mrstats "$harm_dir/rish_l4.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r6=$(mrstats "$harm_dir/rish_l6.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            local r8=$(mrstats "$harm_dir/rish_l8.mif" -mask "$mask" -output mean 2>/dev/null || echo "NA")
            echo "$subj,harmonized,$r0,$r2,$r4,$r6,$r8" >> "$report_file"
        fi
    done
    
    echo "  Summary saved: $report_file"
    cat "$report_file"
    
    echo ""
    echo "  ✓ QC report complete"
}

#=============================================================================
# MAIN
#=============================================================================

main() {
    echo ""
    echo "============================================================"
    echo "  MRTRIX-RISH: Preprocessing + Harmonization Pipeline"
    echo "============================================================"
    echo ""
    
    load_config
    show_config
    
    # Validate paths
    if [[ ! -d "$INPUT_PATH" ]]; then
        echo "ERROR: Input path not found: $INPUT_PATH"
        exit 1
    fi
    
    mkdir -p "$DER_PATH"
    
    # Get all subjects
    local all_subjects="$REFERENCE_SUBJECTS $TARGET_SUBJECTS"
    
    # Determine steps to run
    local run_preprocess=false
    local run_fod=false
    local run_rish=false
    local run_qc=false
    
    if [[ "$STEPS" == "all" ]]; then
        run_preprocess=true
        run_fod=true
        run_rish=true
        run_qc=true
    else
        [[ "$STEPS" == *"preprocess"* ]] && run_preprocess=true
        [[ "$STEPS" == *"fod"* ]] && run_fod=true
        [[ "$STEPS" == *"rish"* ]] && run_rish=true
        [[ "$STEPS" == *"qc"* ]] && run_qc=true
    fi
    
    # Run preprocessing
    if $run_preprocess; then
        for subject in $all_subjects; do
            preprocess_subject "$subject"
        done
    fi
    
    # Run FOD computation
    if $run_fod; then
        for subject in $all_subjects; do
            compute_fod_subject "$subject"
        done
    fi
    
    # Run RISH extraction and harmonization
    if $run_rish; then
        # Extract RISH for all subjects
        for subject in $all_subjects; do
            extract_rish "$subject"
        done
        
        # Create template from reference subjects
        create_rish_template
        
        # Harmonize target subjects
        for subject in $TARGET_SUBJECTS; do
            harmonize_subject "$subject"
        done
    fi
    
    # Generate QC report
    if $run_qc; then
        generate_qc_report
    fi
    
    echo ""
    echo "============================================================"
    echo "  ✓ PIPELINE COMPLETE"
    echo "============================================================"
    echo ""
    echo "Output directory: $DER_PATH"
    echo ""
}

# Run if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi

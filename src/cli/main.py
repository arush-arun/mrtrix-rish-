#!/usr/bin/env python3
"""
mrtrix-rish CLI

MRtrix3-native RISH harmonization pipeline.

Usage:
    mrtrix-rish create-template ...
    mrtrix-rish harmonize ...
    mrtrix-rish qc ...
"""

import argparse
import sys
from pathlib import Path


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser."""
    parser = argparse.ArgumentParser(
        prog="mrtrix-rish",
        description="MRtrix3-native RISH harmonization for multi-site dMRI",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--version", "-V",
        action="version",
        version="%(prog)s 0.1.0"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # create-template subcommand
    template_parser = subparsers.add_parser(
        "create-template",
        help="Create RISH template from reference site"
    )
    template_parser.add_argument(
        "--reference-list", "-r",
        required=True,
        help="Text file listing reference SH images (one per line)"
    )
    template_parser.add_argument(
        "--mask-list", "-m",
        help="Text file listing brain masks (one per line)"
    )
    template_parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for template"
    )
    template_parser.add_argument(
        "--lmax", "-l",
        type=int,
        default=8,
        help="Maximum SH order (default: 8)"
    )
    template_parser.add_argument(
        "--nthreads", "-n",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    template_parser.add_argument(
        "--participants",
        help="Participants TSV/CSV with subject demographics"
    )
    template_parser.add_argument(
        "--covariates",
        help="Comma-separated covariate column names (e.g., age,sex)"
    )
    template_parser.add_argument(
        "--subject-column",
        default="participant_id",
        help="Subject ID column name in participants file (default: participant_id)"
    )

    # harmonize subcommand
    harm_parser = subparsers.add_parser(
        "harmonize",
        help="Harmonize target data using template"
    )
    harm_parser.add_argument(
        "--target", "-t",
        required=True,
        help="Target SH image to harmonize"
    )
    harm_parser.add_argument(
        "--mask", "-m",
        help="Brain mask for target"
    )
    harm_parser.add_argument(
        "--template", "-T",
        required=True,
        help="Template directory (from create-template)"
    )
    harm_parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory"
    )
    harm_parser.add_argument(
        "--config", "-c",
        help="YAML configuration file"
    )
    harm_parser.add_argument(
        "--smoothing-fwhm",
        type=float,
        default=3.0,
        help="Scale map smoothing FWHM in mm (default: 3.0)"
    )
    harm_parser.add_argument(
        "--nthreads", "-n",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    harm_parser.add_argument(
        "--participants",
        help="Participants TSV/CSV for target subject covariates"
    )
    harm_parser.add_argument(
        "--subject-id",
        help="Subject ID to look up in participants file"
    )
    harm_parser.add_argument(
        "--subject-covariates",
        help='Target subject covariates as key=value pairs (e.g., "age=35,sex=M")'
    )

    # qc subcommand
    qc_parser = subparsers.add_parser(
        "qc",
        help="Generate QC report"
    )
    qc_parser.add_argument(
        "--original", "-i",
        required=True,
        help="Original (pre-harmonization) image"
    )
    qc_parser.add_argument(
        "--harmonized", "-H",
        required=True,
        help="Harmonized image"
    )
    qc_parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for QC report"
    )
    qc_parser.add_argument(
        "--mask", "-m",
        help="Brain mask"
    )
    
    # extract-rish subcommand (utility)
    rish_parser = subparsers.add_parser(
        "extract-rish",
        help="Extract RISH features from SH image"
    )
    rish_parser.add_argument(
        "input",
        help="Input SH image"
    )
    rish_parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory"
    )
    rish_parser.add_argument(
        "--mask", "-m",
        help="Brain mask"
    )
    rish_parser.add_argument(
        "--lmax", "-l",
        type=int,
        help="Maximum SH order (auto-detected if not specified)"
    )
    rish_parser.add_argument(
        "--nthreads", "-n",
        type=int,
        default=1,
        help="Number of threads (default: 1)"
    )
    
    # compute-fod subcommand
    fod_parser = subparsers.add_parser(
        "compute-fod",
        help="Compute FOD (auto-detects single/multi-shell)"
    )
    fod_parser.add_argument(
        "dwi",
        help="Input DWI image"
    )
    fod_parser.add_argument(
        "--mask", "-m",
        required=True,
        help="Brain mask"
    )
    fod_parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory"
    )
    fod_parser.add_argument(
        "--algorithm", "-a",
        choices=["auto", "csd", "msmt_csd"],
        default="auto",
        help="FOD algorithm (default: auto-detect)"
    )
    fod_parser.add_argument(
        "--lmax", "-l",
        type=int,
        help="Maximum SH order for FOD"
    )
    fod_parser.add_argument(
        "--nthreads", "-n",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    
    # detect-shells subcommand (utility)
    shells_parser = subparsers.add_parser(
        "detect-shells",
        help="Detect shell structure of DWI data"
    )
    shells_parser.add_argument(
        "dwi",
        help="Input DWI image"
    )
    
    # bids-process subcommand
    bids_parser = subparsers.add_parser(
        "bids",
        help="Process a BIDS dataset"
    )
    bids_parser.add_argument(
        "bids_dir",
        help="Path to BIDS dataset"
    )
    bids_parser.add_argument(
        "--output", "-o",
        default="derivatives",
        help="Output directory (default: BIDS derivatives)"
    )
    bids_parser.add_argument(
        "--subjects", "-s",
        nargs="+",
        help="Specific subjects to process (e.g., sub-01 sub-02)"
    )
    bids_parser.add_argument(
        "--skip-fod",
        action="store_true",
        help="Skip FOD computation (assume input is already SH/FOD)"
    )
    bids_parser.add_argument(
        "--lmax", "-l",
        type=int,
        default=8,
        help="Maximum SH order (default: 8)"
    )
    bids_parser.add_argument(
        "--nthreads", "-n",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    
    # bids-list subcommand
    bids_list_parser = subparsers.add_parser(
        "bids-list",
        help="List DWI files in a BIDS dataset"
    )
    bids_list_parser.add_argument(
        "bids_dir",
        help="Path to BIDS dataset"
    )

    # site-effect subcommand
    site_effect_parser = subparsers.add_parser(
        "site-effect",
        help="Test for site effects using GLM with permutation inference"
    )
    site_effect_parser.add_argument(
        "--site-list", "-s",
        required=True,
        help="CSV file with columns: subject,site,image_path[,age,sex,...]"
    )
    site_effect_parser.add_argument(
        "--mask", "-m",
        required=True,
        help="Brain mask in template space"
    )
    site_effect_parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory for results"
    )
    site_effect_parser.add_argument(
        "--test", "-t",
        choices=["parametric", "permutation"],
        default="permutation",
        help="Test type (default: permutation)"
    )
    site_effect_parser.add_argument(
        "--n-permutations", "-n",
        type=int,
        default=5000,
        help="Number of permutations (default: 5000)"
    )
    site_effect_parser.add_argument(
        "--covariates", "-c",
        help="Comma-separated covariate column names (e.g., age,sex)"
    )
    site_effect_parser.add_argument(
        "--alpha", "-a",
        type=float,
        default=0.05,
        help="Significance level (default: 0.05)"
    )
    site_effect_parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for reproducibility"
    )
    site_effect_parser.add_argument(
        "--heteroscedastic",
        action="store_true",
        help="Use heteroscedastic (G-statistic) test"
    )

    # rish-glm subcommand
    glm_parser = subparsers.add_parser(
        "rish-glm",
        help="RISH-GLM harmonization (joint GLM for site + covariates)"
    )
    glm_parser.add_argument(
        "--manifest", "-i",
        required=True,
        help="CSV with columns: subject,site,fod_path[,covariate columns...]"
    )
    glm_parser.add_argument(
        "--reference-site", "-r",
        required=True,
        help="Name of the reference site"
    )
    glm_parser.add_argument(
        "--mask", "-m",
        required=True,
        help="Brain mask in template space"
    )
    glm_parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory"
    )
    glm_parser.add_argument(
        "--lmax", "-l",
        type=int,
        default=8,
        help="Maximum SH order (default: 8)"
    )
    glm_parser.add_argument(
        "--covariates", "-c",
        help="Comma-separated covariate column names (e.g., age,sex)"
    )
    glm_parser.add_argument(
        "--smoothing-fwhm",
        type=float,
        default=3.0,
        help="Scale map smoothing FWHM in mm (default: 3.0)"
    )
    glm_parser.add_argument(
        "--clip-range",
        default="0.5,2.0",
        help="Scale factor clipping range (default: 0.5,2.0)"
    )
    glm_parser.add_argument(
        "--nthreads", "-n",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    glm_parser.add_argument(
        "--harmonize",
        action="store_true",
        help="Apply harmonization to all target site subjects"
    )

    return parser


def cmd_create_template(args):
    """Handle create-template command."""
    from ..core.rish_features import extract_rish_features
    from ..core.harmonize import run_mrtrix_cmd

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read subject lists
    with open(args.reference_list) as f:
        sh_images = [line.strip() for line in f if line.strip()]

    masks = None
    if args.mask_list:
        with open(args.mask_list) as f:
            masks = [line.strip() for line in f if line.strip()]

    print(f"Creating template from {len(sh_images)} subjects...")

    # Extract RISH for each subject
    rish_list = []
    for i, sh_img in enumerate(sh_images):
        print(f"  Processing subject {i+1}/{len(sh_images)}...")
        mask = masks[i] if masks else None
        subj_dir = output_dir / "subjects" / f"sub-{i:03d}"

        rish = extract_rish_features(
            sh_img,
            str(subj_dir),
            lmax=args.lmax,
            mask=mask,
            n_threads=args.nthreads
        )
        rish_list.append(rish)

    # Covariate adjustment (if participants + covariates provided)
    rish_to_average = rish_list
    if args.participants and args.covariates:
        from ..io.participants import load_participants_tsv, load_participants_csv
        from ..core.covariates import fit_covariate_model, adjust_rish_features_batch

        cov_names = [c.strip() for c in args.covariates.split(",")]
        subj_col = args.subject_column

        # Auto-detect TSV vs CSV
        if args.participants.endswith(".tsv"):
            participants = load_participants_tsv(
                args.participants, cov_names, subject_column=subj_col
            )
        else:
            participants = load_participants_csv(
                args.participants, cov_names, subject_column=subj_col
            )

        if participants.n_subjects != len(sh_images):
            print(
                f"Warning: {participants.n_subjects} participants vs "
                f"{len(sh_images)} images. Ensure ordering matches."
            )

        print(f"  Fitting covariate model ({', '.join(cov_names)})...")

        # Reorganize rish_list into {order: [sub0_path, sub1_path, ...]}
        orders = sorted(rish_list[0].keys())
        rish_by_order = {
            l: [rish_list[i][l] for i in range(len(rish_list))]
            for l in orders
        }

        # Use first mask for fitting (or None)
        fit_mask = masks[0] if masks else None

        model = fit_covariate_model(
            rish_by_order, participants.covariates, fit_mask, str(output_dir)
        )

        print(f"  Adjusting RISH features for {len(sh_images)} subjects...")
        adjusted_list = adjust_rish_features_batch(
            rish_list, participants.covariates, model,
            str(output_dir / "adjusted"), mask_path=fit_mask
        )
        rish_to_average = adjusted_list

    # Average RISH features
    orders = list(rish_to_average[0].keys())
    print(f"Averaging RISH features across subjects...")

    for l in orders:
        images = [r[l] for r in rish_to_average]
        avg_output = output_dir / f"template_rish_l{l}.mif"

        run_mrtrix_cmd([
            "mrmath", *images,
            "mean", str(avg_output),
            "-force",
            "-nthreads", str(args.nthreads)
        ])
        print(f"  Created {avg_output}")

    print(f"✓ Template created in {output_dir}")


def cmd_harmonize(args):
    """Handle harmonize command."""
    from ..core.rish_features import extract_rish_features
    from ..core.scale_maps import compute_scale_maps
    from ..core.harmonize import harmonize_sh
    from ..io.config_io import load_config

    # Load config if provided
    config = None
    if args.config:
        config = load_config(args.config)

    template_dir = Path(args.template)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load template RISH features
    print("Loading template RISH features...")
    template_rish = {}
    for rish_file in template_dir.glob("template_rish_l*.mif"):
        l = int(rish_file.stem.split("_l")[1])
        template_rish[l] = str(rish_file)

    if not template_rish:
        print(f"Error: No template RISH files found in {template_dir}")
        sys.exit(1)

    lmax = max(template_rish.keys())
    print(f"  Found template with lmax={lmax}")

    # Extract target RISH
    print("Extracting target RISH features...")
    target_rish_dir = output_dir / "target_rish"
    target_rish = extract_rish_features(
        args.target,
        str(target_rish_dir),
        lmax=lmax,
        mask=args.mask,
        n_threads=args.nthreads
    )

    # Check for covariate model in template directory
    covariate_model_path = template_dir / "covariate_model" / "covariate_model.json"
    rish_for_scale = target_rish

    if covariate_model_path.exists():
        from ..core.covariates import load_covariate_model, adjust_rish_features

        model = load_covariate_model(str(covariate_model_path))
        print(f"  Loaded covariate model ({', '.join(model.covariate_names)})")

        # Parse target subject covariates
        subject_covs = _parse_target_covariates(args, model.covariate_names)

        if subject_covs:
            print(f"  Adjusting target RISH for covariates...")
            adjusted_dir = output_dir / "target_rish_adjusted"
            adjusted_rish = adjust_rish_features(
                target_rish, subject_covs, model,
                str(adjusted_dir), mask_path=args.mask
            )
            rish_for_scale = adjusted_rish
        else:
            print(
                "  Warning: Covariate model found but no subject covariates "
                "provided. Use --subject-covariates or --participants + "
                "--subject-id. Proceeding without adjustment."
            )

    # Compute scale maps
    print("Computing scale maps...")
    scale_dir = output_dir / "scale_maps"
    scale_maps = compute_scale_maps(
        template_rish,
        rish_for_scale,
        str(scale_dir),
        mask=args.mask,
        smoothing_fwhm=args.smoothing_fwhm,
        n_threads=args.nthreads
    )

    # Apply harmonization to ORIGINAL target SH (not adjusted)
    print("Applying harmonization...")
    harmonized = output_dir / "sh_harmonized.mif"
    harmonize_sh(
        args.target,
        scale_maps,
        str(harmonized),
        lmax=lmax,
        n_threads=args.nthreads
    )

    print(f"✓ Harmonized output: {harmonized}")


def _parse_target_covariates(args, covariate_names):
    """Parse target subject covariates from CLI args.

    Supports two modes:
    1. --subject-covariates "age=35,sex=M"
    2. --participants file.tsv --subject-id sub-01

    Returns
    -------
    dict or None
        covariate_name -> raw value, or None if not provided
    """
    from ..io.participants import _encode_categorical

    if args.subject_covariates:
        # Parse "age=35,sex=M" format
        covs = {}
        for pair in args.subject_covariates.split(","):
            key, val = pair.strip().split("=", 1)
            key = key.strip()
            val = val.strip()

            # Try numeric first
            try:
                covs[key] = float(val)
            except ValueError:
                # Categorical: encode sex M/F
                covs[key] = _encode_categorical([val])[0]

        return covs

    if args.participants and args.subject_id:
        from ..io.participants import load_participants_tsv, load_participants_csv

        if args.participants.endswith(".tsv"):
            data = load_participants_tsv(
                args.participants, covariate_names,
                subject_ids=[args.subject_id]
            )
        else:
            data = load_participants_csv(
                args.participants, covariate_names,
                subject_ids=[args.subject_id]
            )

        return {name: vals[0] for name, vals in data.covariates.items()}

    return None


def cmd_extract_rish(args):
    """Handle extract-rish command."""
    from ..core.rish_features import extract_rish_features
    
    print(f"Extracting RISH features from {args.input}...")
    
    rish = extract_rish_features(
        args.input,
        args.output,
        lmax=args.lmax,
        mask=args.mask,
        n_threads=args.nthreads
    )
    
    print("Output files:")
    for l, path in sorted(rish.items()):
        print(f"  l={l}: {path}")


def cmd_qc(args):
    """Handle qc command."""
    print("QC report generation not yet implemented.")
    print("Coming soon!")


def cmd_compute_fod(args):
    """Handle compute-fod command."""
    from ..core.fod import compute_fod, detect_shells
    
    print(f"Analyzing {args.dwi}...")
    shell_info = detect_shells(args.dwi)
    
    print(f"\nShell structure detected:")
    print(f"  Type: {shell_info.shell_type.value}")
    print(f"  b-values: {shell_info.bvalues}")
    print(f"  Volumes per shell: {shell_info.shell_counts}")
    print(f"  b=0 volumes: {shell_info.n_b0}")
    
    algorithm = args.algorithm
    if algorithm == "auto":
        algorithm = "msmt_csd" if shell_info.is_multi_shell else "csd"
    print(f"\nUsing algorithm: {algorithm}")
    
    result = compute_fod(
        args.dwi,
        args.mask,
        args.output,
        shell_info=shell_info,
        lmax=args.lmax,
        algorithm=args.algorithm,
        n_threads=args.nthreads
    )
    
    print(f"\n✓ FOD computed:")
    print(f"  WM FOD: {result['wm_fod']}")
    if 'gm' in result:
        print(f"  GM: {result['gm']}")
        print(f"  CSF: {result['csf']}")


def cmd_detect_shells(args):
    """Handle detect-shells command."""
    from ..core.fod import detect_shells
    
    shell_info = detect_shells(args.dwi)
    
    print(f"DWI: {args.dwi}")
    print(f"\nShell structure:")
    print(f"  Type: {shell_info.shell_type.value}")
    print(f"  b-values: {shell_info.bvalues}")
    print(f"  Volumes per shell: {shell_info.shell_counts}")
    print(f"  b=0 volumes: {shell_info.n_b0}")
    print(f"  Total volumes: {shell_info.n_total}")
    
    if shell_info.is_multi_shell:
        print(f"\n→ Recommended: MSMT-CSD (dwi2fod msmt_csd)")
    else:
        print(f"\n→ Recommended: CSD (dwi2fod csd)")


def cmd_bids(args):
    """Handle bids command."""
    from ..core.bids_workflow import process_bids_dataset
    
    process_bids_dataset(
        bids_dir=args.bids_dir,
        output_dir=args.output,
        subjects=args.subjects,
        compute_fod_flag=not args.skip_fod,
        lmax=args.lmax,
        n_threads=args.nthreads
    )


def cmd_bids_list(args):
    """Handle bids-list command."""
    from ..io.bids_io import find_bids_dwi

    entries = find_bids_dwi(args.bids_dir)

    if not entries:
        print(f"No DWI files found in {args.bids_dir}")
        return

    print(f"Found {len(entries)} DWI datasets in {args.bids_dir}:\n")

    for e in entries:
        session_str = f" / {e['session']}" if e.get('session') else ""
        print(f"  {e['subject']}{session_str}")
        print(f"    DWI:  {e['dwi']}")
        print(f"    bval: {e['bval'] or '(not found)'}")
        print(f"    bvec: {e['bvec'] or '(not found)'}")
        print(f"    mask: {e['mask'] or '(not found)'}")
        print()


def cmd_rish_glm(args):
    """Handle rish-glm command."""
    import csv
    from ..core.rish_features import extract_rish_features
    from ..core.rish_glm import (
        fit_rish_glm,
        compute_glm_scale_maps,
    )
    from ..core.harmonize import harmonize_sh
    from ..io.participants import _encode_categorical

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse manifest CSV
    print(f"Reading manifest from {args.manifest}...")
    subjects = []
    site_labels = []
    fod_paths = []
    covariates = {}
    covariate_names = []

    if args.covariates:
        covariate_names = [c.strip() for c in args.covariates.split(",")]
        for name in covariate_names:
            covariates[name] = []

    with open(args.manifest, "r") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []

        # Validate required columns
        for col in ["subject", "site", "fod_path"]:
            if col not in fieldnames:
                print(f"Error: manifest must have '{col}' column")
                sys.exit(1)

        for row in reader:
            subjects.append(row["subject"])
            site_labels.append(row["site"])
            fod_paths.append(row["fod_path"])

            for name in covariate_names:
                if name not in row:
                    print(f"Error: covariate column '{name}' not found in manifest")
                    sys.exit(1)
                val = row[name]
                try:
                    covariates[name].append(float(val))
                except ValueError:
                    covariates[name].append(
                        _encode_categorical([val])[0]
                    )

    unique_sites = sorted(set(site_labels))
    n_subjects = len(subjects)

    print(f"\nData summary:")
    for site in unique_sites:
        n = sum(1 for s in site_labels if s == site)
        print(f"  {site}: {n} subjects")
    if covariate_names:
        print(f"  Covariates: {', '.join(covariate_names)}")
    print(f"  Reference site: {args.reference_site}")

    # Extract RISH features for all subjects
    print(f"\nExtracting RISH features...")
    all_rish = []
    for i, fod_path in enumerate(fod_paths):
        subj_dir = output_dir / "rish" / subjects[i]
        print(f"  [{i+1}/{n_subjects}] {subjects[i]}...")
        rish = extract_rish_features(
            fod_path,
            str(subj_dir),
            lmax=args.lmax,
            mask=args.mask,
            n_threads=args.nthreads,
        )
        all_rish.append(rish)

    # Organize RISH by order
    orders = sorted(all_rish[0].keys())
    rish_by_order = {
        l: [all_rish[i][l] for i in range(n_subjects)]
        for l in orders
    }

    # Fit RISH-GLM
    print(f"\nFitting RISH-GLM...")
    result = fit_rish_glm(
        rish_by_order,
        site_labels,
        mask_path=args.mask,
        output_dir=str(output_dir),
        reference_site=args.reference_site,
        covariates=covariates if covariates else None,
    )

    print(f"  Model fitted with {len(result.design_columns)} columns: "
          f"{', '.join(result.design_columns)}")

    # Compute scale maps for each target site
    clip_range = tuple(float(x) for x in args.clip_range.split(","))
    target_sites = [s for s in unique_sites if s != args.reference_site]

    for target in target_sites:
        print(f"\nComputing scale maps for {target}...")
        scale_dir = output_dir / "scale_maps" / target
        scale_maps = compute_glm_scale_maps(
            result,
            target,
            str(scale_dir),
            mask=args.mask,
            smoothing_fwhm=args.smoothing_fwhm,
            clip_range=clip_range,
            n_threads=args.nthreads,
        )
        print(f"  Created {len(scale_maps)} scale maps")

        # Apply harmonization if requested
        if args.harmonize:
            print(f"  Harmonizing {target} subjects...")
            harm_dir = output_dir / "harmonized" / target
            harm_dir.mkdir(parents=True, exist_ok=True)

            for i, (subj, site) in enumerate(zip(subjects, site_labels)):
                if site != target:
                    continue
                out_path = harm_dir / f"{subj}_harmonized.mif"
                harmonize_sh(
                    fod_paths[i],
                    scale_maps,
                    str(out_path),
                    lmax=args.lmax,
                    n_threads=args.nthreads,
                )
                print(f"    {subj} -> {out_path}")

    print(f"\n✓ RISH-GLM complete. Output: {output_dir}")


def cmd_site_effect(args):
    """Handle site-effect command."""
    import csv
    from ..qc.site_effects import test_site_effect

    # Parse CSV file
    print(f"Reading site list from {args.site_list}...")
    image_paths = {}
    covariates = {}
    covariate_names = []

    if args.covariates:
        covariate_names = [c.strip() for c in args.covariates.split(",")]
        for name in covariate_names:
            covariates[name] = []

    with open(args.site_list, "r") as f:
        reader = csv.DictReader(f)

        # Validate required columns
        fieldnames = reader.fieldnames or []
        if "subject" not in fieldnames or "site" not in fieldnames:
            print("Error: CSV must have 'subject' and 'site' columns")
            sys.exit(1)

        # Check for image path column
        image_col = None
        for col in ["image_path", "image", "path", "fa_path", "fa"]:
            if col in fieldnames:
                image_col = col
                break
        if image_col is None:
            print("Error: CSV must have an image path column (image_path, image, path, fa_path, or fa)")
            sys.exit(1)

        # Read data
        for row in reader:
            site = row["site"]
            img_path = row[image_col]

            if site not in image_paths:
                image_paths[site] = []
            image_paths[site].append(img_path)

            # Read covariates
            for name in covariate_names:
                if name in row:
                    val = row[name]
                    # Handle sex as numeric
                    if name.lower() == "sex":
                        val = 1.0 if val.upper() in ["M", "MALE", "1"] else 0.0
                    else:
                        val = float(val)
                    covariates[name].append(val)

    # Print summary
    print(f"\nData summary:")
    for site, paths in sorted(image_paths.items()):
        print(f"  {site}: {len(paths)} subjects")
    if covariate_names:
        print(f"  Covariates: {', '.join(covariate_names)}")

    # Set up variance groups if heteroscedastic
    variance_groups = None
    if args.heteroscedastic:
        # Map sites to variance group indices
        unique_sites = sorted(image_paths.keys())
        variance_groups = {site: i for i, site in enumerate(unique_sites)}
        print(f"  Using heteroscedastic test with {len(unique_sites)} variance groups")

    # Run site effect test
    print(f"\nRunning site effect analysis...")
    result = test_site_effect(
        image_paths=image_paths,
        mask_path=args.mask,
        output_dir=args.output,
        covariates=covariates if covariates else None,
        n_permutations=args.n_permutations,
        alpha=args.alpha,
        variance_groups=variance_groups,
        seed=args.seed,
        save_maps=True,
        verbose=True
    )

    # Print conclusion
    print("\n" + "=" * 60)
    if result.percent_significant_permutation < 5.0:
        print("CONCLUSION: No significant site effect detected")
        print("            (Harmonization successful)")
    elif result.percent_significant_permutation < 15.0:
        print("CONCLUSION: Moderate site effect detected")
        print("            (Harmonization partially successful)")
    else:
        print("CONCLUSION: Significant site effect detected")
        print("            (Harmonization may be needed)")
    print("=" * 60)


def main():
    """Main entry point."""
    parser = create_parser()
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    
    # Dispatch to command handler
    handlers = {
        "create-template": cmd_create_template,
        "harmonize": cmd_harmonize,
        "extract-rish": cmd_extract_rish,
        "compute-fod": cmd_compute_fod,
        "detect-shells": cmd_detect_shells,
        "bids": cmd_bids,
        "bids-list": cmd_bids_list,
        "qc": cmd_qc,
        "site-effect": cmd_site_effect,
        "rish-glm": cmd_rish_glm,
    }
    
    handler = handlers.get(args.command)
    if handler:
        handler(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()

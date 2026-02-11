"""
BIDS Workflow Module

High-level workflows for processing BIDS datasets.
"""

from pathlib import Path
from typing import List, Dict, Optional
import tempfile
import json

from .fod import compute_fod, detect_shells
from .rish_features import extract_rish_features
from .scale_maps import compute_scale_maps
from .harmonize import harmonize_sh
from ..io.bids_io import find_bids_dwi, convert_to_mif, get_bids_derivative_path


def process_bids_subject(
    dwi_entry: Dict,
    output_root: str,
    compute_fod_flag: bool = True,
    lmax: int = 8,
    n_threads: int = 4
) -> Dict[str, str]:
    """
    Process a single BIDS subject through FOD and RISH extraction.
    
    Parameters
    ----------
    dwi_entry : dict
        Entry from find_bids_dwi() with dwi, bval, bvec, mask paths
    output_root : str
        Output directory root
    compute_fod_flag : bool
        Whether to compute FOD (if False, assumes input is already SH/FOD)
    lmax : int
        Maximum SH order
    n_threads : int
        Number of threads
        
    Returns
    -------
    dict
        Paths to outputs (mif, fod, rish, etc.)
    """
    subject = dwi_entry["subject"]
    session = dwi_entry.get("session")
    
    output_dir = Path(output_root)
    if session:
        subj_dir = output_dir / subject / session
    else:
        subj_dir = output_dir / subject
    subj_dir.mkdir(parents=True, exist_ok=True)
    
    outputs = {"subject": subject, "session": session}
    
    # Convert to MIF if needed
    dwi_path = dwi_entry["dwi"]
    if dwi_path.endswith((".nii", ".nii.gz")):
        mif_path = subj_dir / "dwi.mif"
        if not mif_path.exists():
            print(f"  Converting {subject} to MIF format...")
            convert_to_mif(
                dwi_path,
                dwi_entry["bval"],
                dwi_entry["bvec"],
                str(mif_path),
                dwi_entry.get("json")
            )
        dwi_path = str(mif_path)
        outputs["mif"] = dwi_path
    
    # Get or create mask
    mask_path = dwi_entry.get("mask")
    if not mask_path:
        mask_path = subj_dir / "mask.mif"
        if not mask_path.exists():
            print(f"  Creating brain mask for {subject}...")
            import subprocess
            subprocess.run([
                "dwi2mask", dwi_path, str(mask_path), "-force"
            ], check=True, capture_output=True)
        mask_path = str(mask_path)
    outputs["mask"] = mask_path
    
    # Compute FOD
    if compute_fod_flag:
        fod_dir = subj_dir / "fod"
        print(f"  Computing FOD for {subject}...")
        fod_result = compute_fod(
            dwi_path,
            mask_path,
            str(fod_dir),
            lmax=lmax,
            n_threads=n_threads
        )
        outputs["fod"] = fod_result
        sh_image = fod_result["wm_fod"]
    else:
        sh_image = dwi_path
    
    outputs["sh"] = sh_image
    
    # Extract RISH features
    rish_dir = subj_dir / "rish"
    print(f"  Extracting RISH features for {subject}...")
    rish = extract_rish_features(
        sh_image,
        str(rish_dir),
        lmax=lmax,
        mask=mask_path,
        n_threads=n_threads
    )
    outputs["rish"] = rish
    
    return outputs


def process_bids_dataset(
    bids_dir: str,
    output_dir: str,
    subjects: Optional[List[str]] = None,
    compute_fod_flag: bool = True,
    lmax: int = 8,
    n_threads: int = 4
) -> List[Dict]:
    """
    Process all subjects in a BIDS dataset.
    
    Parameters
    ----------
    bids_dir : str
        Path to BIDS dataset
    output_dir : str
        Output directory (or "derivatives" to write to BIDS derivatives)
    subjects : list of str, optional
        Specific subjects to process (default: all)
    compute_fod_flag : bool
        Compute FODs from DWI
    lmax : int
        Maximum SH order
    n_threads : int
        Threads per subject
        
    Returns
    -------
    list of dict
        Processing results for each subject
    """
    bids_dir = Path(bids_dir)
    
    # Handle output directory
    if output_dir == "derivatives":
        output_dir = bids_dir / "derivatives" / "mrtrix-rish"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all DWI files
    print(f"Scanning BIDS dataset: {bids_dir}")
    dwi_entries = find_bids_dwi(str(bids_dir))
    
    if subjects:
        dwi_entries = [e for e in dwi_entries if e["subject"] in subjects]
    
    print(f"Found {len(dwi_entries)} DWI datasets")
    
    results = []
    for i, entry in enumerate(dwi_entries):
        print(f"\n[{i+1}/{len(dwi_entries)}] Processing {entry['subject']}" + 
              (f" / {entry['session']}" if entry.get('session') else ""))
        
        try:
            result = process_bids_subject(
                entry,
                str(output_dir),
                compute_fod_flag=compute_fod_flag,
                lmax=lmax,
                n_threads=n_threads
            )
            result["status"] = "success"
        except Exception as e:
            print(f"  ERROR: {e}")
            result = {
                "subject": entry["subject"],
                "session": entry.get("session"),
                "status": "error",
                "error": str(e)
            }
        
        results.append(result)
    
    # Save processing log
    log_file = output_dir / "processing_log.json"
    with open(log_file, "w") as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n✓ Processed {len(results)} subjects")
    print(f"  Log: {log_file}")
    
    return results


def harmonize_bids_sites(
    reference_outputs: List[Dict],
    target_outputs: List[Dict],
    output_dir: str,
    smoothing_fwhm: float = 3.0,
    n_threads: int = 4
) -> Dict:
    """
    Harmonize target site data to reference site.
    
    Parameters
    ----------
    reference_outputs : list of dict
        Outputs from process_bids_subject for reference site
    target_outputs : list of dict
        Outputs from process_bids_subject for target site
    output_dir : str
        Output directory for harmonized data
    smoothing_fwhm : float
        Scale map smoothing
    n_threads : int
        Number of threads
        
    Returns
    -------
    dict
        Harmonization results including template and harmonized subjects
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get lmax from first reference subject
    first_rish = reference_outputs[0]["rish"]
    lmax = max(first_rish.keys())
    
    # Average reference RISH to create template
    print("Creating reference template...")
    template_dir = output_dir / "template"
    template_dir.mkdir(exist_ok=True)
    
    template_rish = {}
    for l in range(0, lmax + 1, 2):
        ref_images = [r["rish"][l] for r in reference_outputs if r.get("rish")]
        template_file = template_dir / f"template_rish_l{l}.mif"
        
        import subprocess
        subprocess.run([
            "mrmath", *ref_images,
            "mean", str(template_file),
            "-force", "-nthreads", str(n_threads)
        ], check=True, capture_output=True)
        
        template_rish[l] = str(template_file)
    
    # Harmonize each target subject
    harmonized_results = []
    for i, target in enumerate(target_outputs):
        if target.get("status") != "success":
            continue
            
        subj = target["subject"]
        print(f"Harmonizing {subj} ({i+1}/{len(target_outputs)})...")
        
        subj_dir = output_dir / "harmonized" / subj
        subj_dir.mkdir(parents=True, exist_ok=True)
        
        # Compute scale maps
        scale_maps = compute_scale_maps(
            template_rish,
            target["rish"],
            str(subj_dir / "scale_maps"),
            mask=target.get("mask"),
            smoothing_fwhm=smoothing_fwhm,
            n_threads=n_threads
        )
        
        # Apply harmonization
        harmonized_sh = subj_dir / "sh_harmonized.mif"
        harmonize_sh(
            target["sh"],
            scale_maps,
            str(harmonized_sh),
            lmax=lmax,
            n_threads=n_threads
        )
        
        harmonized_results.append({
            "subject": subj,
            "harmonized_sh": str(harmonized_sh),
            "scale_maps": scale_maps
        })
    
    return {
        "template": template_rish,
        "harmonized": harmonized_results,
        "n_reference": len(reference_outputs),
        "n_target": len(harmonized_results)
    }

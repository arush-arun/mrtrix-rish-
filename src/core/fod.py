"""
FOD Computation Module

Automatically detects single-shell vs multi-shell data and computes
FODs using the appropriate algorithm:
- Single-shell: CSD (Constrained Spherical Deconvolution)
- Multi-shell: MSMT-CSD (Multi-Shell Multi-Tissue CSD)
"""

import subprocess
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from enum import Enum
import tempfile


class ShellType(Enum):
    """DWI acquisition type."""
    SINGLE_SHELL = "single_shell"
    MULTI_SHELL = "multi_shell"


@dataclass
class ShellInfo:
    """Information about DWI shells."""
    shell_type: ShellType
    bvalues: List[int]          # Unique b-values (excluding b=0)
    shell_counts: Dict[int, int]  # b-value -> number of volumes
    n_b0: int                    # Number of b=0 volumes
    n_total: int                 # Total volumes
    
    @property
    def n_shells(self) -> int:
        """Number of non-zero shells."""
        return len(self.bvalues)
    
    @property
    def is_multi_shell(self) -> bool:
        return self.shell_type == ShellType.MULTI_SHELL


def run_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command."""
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def detect_shells(
    dwi: str,
    b0_threshold: float = 50.0,
    shell_tolerance: float = 100.0
) -> ShellInfo:
    """
    Detect shell structure from DWI image.
    
    Parameters
    ----------
    dwi : str
        Path to DWI image (must have gradient info in header or sidecar)
    b0_threshold : float
        Maximum b-value to consider as b=0 (default: 50)
    shell_tolerance : float
        Tolerance for grouping b-values into shells (default: 100)
        
    Returns
    -------
    ShellInfo
        Shell structure information
    """
    # Use mrinfo to get b-values
    result = run_cmd(["mrinfo", "-shell_bvalues", dwi])
    bval_str = result.stdout.strip()
    
    if not bval_str:
        # Try to get from -bvalue_scaling
        result = run_cmd(["mrinfo", "-dwgrad", dwi])
        # Parse gradient table... for now raise error
        raise ValueError(f"Could not detect b-values from {dwi}")
    
    # Parse b-values (space-separated)
    all_bvals = [float(b) for b in bval_str.split()]
    
    # Get shell counts
    result = run_cmd(["mrinfo", "-shell_sizes", dwi])
    shell_sizes = [int(s) for s in result.stdout.strip().split()]
    
    # Separate b=0 and non-zero shells
    bvalues = []
    shell_counts = {}
    n_b0 = 0
    
    for bval, count in zip(all_bvals, shell_sizes):
        if bval <= b0_threshold:
            n_b0 = count
        else:
            bval_int = int(round(bval))
            bvalues.append(bval_int)
            shell_counts[bval_int] = count
    
    # Determine shell type
    shell_type = ShellType.MULTI_SHELL if len(bvalues) > 1 else ShellType.SINGLE_SHELL
    
    n_total = n_b0 + sum(shell_counts.values())
    
    return ShellInfo(
        shell_type=shell_type,
        bvalues=sorted(bvalues),
        shell_counts=shell_counts,
        n_b0=n_b0,
        n_total=n_total
    )


def estimate_response(
    dwi: str,
    mask: str,
    output_dir: str,
    shell_info: Optional[ShellInfo] = None,
    algorithm: str = "dhollander",
    n_threads: int = 1
) -> Dict[str, str]:
    """
    Estimate tissue response functions.
    
    Parameters
    ----------
    dwi : str
        DWI image
    mask : str
        Brain mask
    output_dir : str
        Output directory for response files
    shell_info : ShellInfo, optional
        Pre-computed shell info (auto-detected if None)
    algorithm : str
        Response estimation algorithm (dhollander, tournier, tax, fa)
    n_threads : int
        Number of threads
        
    Returns
    -------
    dict
        Paths to response files (wm, gm, csf depending on shell type)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if shell_info is None:
        shell_info = detect_shells(dwi)
    
    thread_opt = ["-nthreads", str(n_threads)]
    
    responses = {}
    
    if shell_info.is_multi_shell:
        # Multi-shell: estimate WM, GM, CSF responses
        wm_response = output_dir / "response_wm.txt"
        gm_response = output_dir / "response_gm.txt"
        csf_response = output_dir / "response_csf.txt"
        
        cmd = [
            "dwi2response", algorithm,
            dwi,
            str(wm_response),
            str(gm_response),
            str(csf_response),
            "-mask", mask,
            "-force"
        ] + thread_opt
        
        # dhollander doesn't need explicit shells, others might
        if algorithm == "dhollander":
            voxels = output_dir / "response_voxels.mif"
            cmd.extend(["-voxels", str(voxels)])
        
        run_cmd(cmd)
        
        responses = {
            "wm": str(wm_response),
            "gm": str(gm_response),
            "csf": str(csf_response)
        }
        
    else:
        # Single-shell: estimate WM response only
        wm_response = output_dir / "response_wm.txt"
        
        # Use tournier algorithm for single-shell
        ss_algorithm = "tournier" if algorithm == "dhollander" else algorithm
        
        cmd = [
            "dwi2response", ss_algorithm,
            dwi,
            str(wm_response),
            "-mask", mask,
            "-force"
        ] + thread_opt
        
        run_cmd(cmd)
        
        responses = {
            "wm": str(wm_response)
        }
    
    return responses


def compute_fod(
    dwi: str,
    mask: str,
    output_dir: str,
    response: Optional[Dict[str, str]] = None,
    shell_info: Optional[ShellInfo] = None,
    lmax: Optional[int] = None,
    algorithm: str = "auto",
    n_threads: int = 1
) -> Dict[str, str]:
    """
    Compute FODs from DWI data.
    
    Automatically selects CSD or MSMT-CSD based on shell structure.
    
    Parameters
    ----------
    dwi : str
        DWI image
    mask : str
        Brain mask
    output_dir : str
        Output directory
    response : dict, optional
        Pre-computed response functions. If None, will estimate.
    shell_info : ShellInfo, optional
        Pre-computed shell info
    lmax : int, optional
        Maximum SH order for FOD (auto if None)
    algorithm : str
        "auto", "csd", or "msmt_csd"
    n_threads : int
        Number of threads
        
    Returns
    -------
    dict
        Paths to FOD images (wm_fod, gm, csf for multi-shell)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    thread_opt = ["-nthreads", str(n_threads)]
    
    # Detect shells if needed
    if shell_info is None:
        shell_info = detect_shells(dwi)
    
    print(f"Detected: {shell_info.shell_type.value}")
    print(f"  Shells: {shell_info.bvalues}")
    print(f"  b=0 volumes: {shell_info.n_b0}")
    
    # Determine algorithm
    if algorithm == "auto":
        algorithm = "msmt_csd" if shell_info.is_multi_shell else "csd"
    
    # Estimate response if not provided
    if response is None:
        response_dir = output_dir / "response"
        response = estimate_response(
            dwi, mask, str(response_dir),
            shell_info=shell_info,
            n_threads=n_threads
        )
    
    outputs = {}
    
    if algorithm == "msmt_csd" and shell_info.is_multi_shell:
        # Multi-Shell Multi-Tissue CSD
        wm_fod = output_dir / "wm_fod.mif"
        gm = output_dir / "gm.mif"
        csf = output_dir / "csf.mif"
        
        cmd = [
            "dwi2fod", "msmt_csd",
            dwi,
            response["wm"], str(wm_fod),
            response["gm"], str(gm),
            response["csf"], str(csf),
            "-mask", mask,
            "-force"
        ] + thread_opt
        
        if lmax:
            cmd.extend(["-lmax", f"{lmax},{lmax},{lmax}"])
        
        run_cmd(cmd)
        
        outputs = {
            "wm_fod": str(wm_fod),
            "gm": str(gm),
            "csf": str(csf),
            "algorithm": "msmt_csd"
        }
        
        # Also create combined tissue image for visualization
        vf = output_dir / "tissue_vf.mif"
        run_cmd([
            "mrcat", str(csf), str(gm), str(wm_fod),
            str(vf), "-axis", "3", "-force"
        ] + thread_opt)
        outputs["tissue_vf"] = str(vf)
        
    else:
        # Standard single-tissue CSD
        wm_fod = output_dir / "wm_fod.mif"
        
        cmd = [
            "dwi2fod", "csd",
            dwi,
            response["wm"], str(wm_fod),
            "-mask", mask,
            "-force"
        ] + thread_opt
        
        if lmax:
            cmd.extend(["-lmax", str(lmax)])
        
        run_cmd(cmd)
        
        outputs = {
            "wm_fod": str(wm_fod),
            "algorithm": "csd"
        }
    
    # Save shell info
    info_file = output_dir / "shell_info.txt"
    with open(info_file, "w") as f:
        f.write(f"shell_type: {shell_info.shell_type.value}\n")
        f.write(f"bvalues: {shell_info.bvalues}\n")
        f.write(f"shell_counts: {shell_info.shell_counts}\n")
        f.write(f"n_b0: {shell_info.n_b0}\n")
        f.write(f"algorithm: {outputs['algorithm']}\n")
    
    outputs["shell_info"] = str(info_file)
    outputs["response"] = response
    
    return outputs


def compute_fod_batch(
    dwi_list: List[str],
    mask_list: List[str],
    output_dir: str,
    n_threads: int = 4
) -> List[Dict[str, str]]:
    """
    Compute FODs for multiple subjects.
    
    Parameters
    ----------
    dwi_list : list of str
        DWI images
    mask_list : list of str
        Brain masks
    output_dir : str
        Base output directory
    n_threads : int
        Threads per subject
        
    Returns
    -------
    list of dict
        FOD outputs for each subject
    """
    output_dir = Path(output_dir)
    results = []
    
    for i, (dwi, mask) in enumerate(zip(dwi_list, mask_list)):
        print(f"Processing subject {i+1}/{len(dwi_list)}...")
        subj_dir = output_dir / f"sub-{i:03d}"
        
        fod = compute_fod(
            dwi, mask, str(subj_dir),
            n_threads=n_threads
        )
        results.append(fod)
    
    return results


def fod_to_sh(fod: str, output: str, lmax: int = 8) -> str:
    """
    Extract SH coefficients from FOD for RISH harmonization.
    
    For WM FOD, this is typically already in SH format.
    This function ensures correct lmax and format.
    
    Parameters
    ----------
    fod : str
        FOD image
    output : str
        Output SH image
    lmax : int
        Target lmax
        
    Returns
    -------
    str
        Path to SH image
    """
    # Check current lmax
    result = run_cmd(["mrinfo", "-size", fod])
    sizes = result.stdout.strip().split()
    n_vols = int(sizes[3]) if len(sizes) >= 4 else int(sizes[-1])
    
    # Calculate current lmax
    import math
    current_lmax = int((-3 + math.sqrt(1 + 8 * n_vols)) / 2)
    
    if current_lmax == lmax:
        # Just copy
        run_cmd(["mrconvert", fod, output, "-force"])
    elif current_lmax > lmax:
        # Truncate to lower lmax
        n_target = (lmax + 1) * (lmax + 2) // 2
        run_cmd([
            "mrconvert", fod,
            "-coord", "3", f"0:{n_target}",
            output, "-force"
        ])
    else:
        # Pad with zeros (unusual case)
        print(f"Warning: FOD lmax ({current_lmax}) < target ({lmax}), padding with zeros")
        run_cmd(["mrconvert", fod, output, "-force"])
    
    return output


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python fod.py <dwi> <mask> [output_dir]")
        sys.exit(1)
    
    dwi = sys.argv[1]
    mask = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else "fod_output"
    
    print(f"Computing FOD for {dwi}...")
    result = compute_fod(dwi, mask, output_dir)
    print(f"Output: {result}")

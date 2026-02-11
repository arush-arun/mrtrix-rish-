"""
Template Creation Module

Creates population templates and registers subjects to common space.
"""

import subprocess
from pathlib import Path
from typing import List, Optional
import tempfile


def run_mrtrix_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run an MRtrix3 command."""
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def create_template(
    fod_images: List[str],
    mask_images: List[str],
    output_dir: str,
    initial_template: Optional[str] = None,
    n_threads: int = 4
) -> dict:
    """
    Create a population template from FOD images.
    
    Uses MRtrix3's population_template for FOD-based registration.
    
    Parameters
    ----------
    fod_images : list of str
        List of FOD images from all subjects (all sites)
    mask_images : list of str
        Corresponding brain masks
    output_dir : str
        Output directory
    initial_template : str, optional
        Initial template to register to (e.g., MNI)
    n_threads : int
        Number of threads
        
    Returns
    -------
    dict
        Paths to template outputs
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create input file lists
    fod_list = output_dir / "fod_list.txt"
    mask_list = output_dir / "mask_list.txt"
    
    with open(fod_list, "w") as f:
        f.write("\n".join(fod_images))
    
    with open(mask_list, "w") as f:
        f.write("\n".join(mask_images))
    
    template_fod = output_dir / "template_fod.mif"
    template_mask = output_dir / "template_mask.mif"
    warps_dir = output_dir / "warps"
    warps_dir.mkdir(exist_ok=True)
    
    # Build population_template command
    cmd = [
        "population_template",
        str(fod_list),
        str(template_fod),
        "-mask_dir", str(mask_list),
        "-template_mask", str(template_mask),
        "-warp_dir", str(warps_dir),
        "-nthreads", str(n_threads)
    ]
    
    if initial_template:
        cmd.extend(["-template", initial_template])
    
    run_mrtrix_cmd(cmd)
    
    return {
        "template_fod": str(template_fod),
        "template_mask": str(template_mask),
        "warps_dir": str(warps_dir),
        "fod_list": str(fod_list),
        "mask_list": str(mask_list)
    }


def warp_to_template(
    image: str,
    warp: str,
    output: str,
    template: str,
    interp: str = "linear",
    n_threads: int = 1
) -> str:
    """
    Warp an image to template space.
    
    Parameters
    ----------
    image : str
        Input image to warp
    warp : str
        Warp field (from population_template)
    output : str
        Output warped image
    template : str
        Template image (for grid)
    interp : str
        Interpolation method
    n_threads : int
        Number of threads
        
    Returns
    -------
    str
        Path to warped image
    """
    run_mrtrix_cmd([
        "mrtransform", image,
        "-warp", warp,
        "-template", template,
        "-interp", interp,
        output,
        "-force",
        "-nthreads", str(n_threads)
    ])
    
    return output


def dwi_to_sh(
    dwi: str,
    output: str,
    lmax: int = 8,
    mask: Optional[str] = None,
    n_threads: int = 1
) -> str:
    """
    Convert DWI to SH coefficients.
    
    Parameters
    ----------
    dwi : str
        Input DWI image
    output : str
        Output SH image
    lmax : int
        Maximum SH order
    mask : str, optional
        Brain mask
    n_threads : int
        Number of threads
        
    Returns
    -------
    str
        Path to SH image
    """
    cmd = [
        "amp2sh", dwi,
        "-lmax", str(lmax),
        output,
        "-force",
        "-nthreads", str(n_threads)
    ]
    
    if mask:
        cmd.extend(["-mask", mask])
    
    run_mrtrix_cmd(cmd)
    
    return output


def sh_to_dwi(
    sh_image: str,
    directions: str,
    output: str,
    n_threads: int = 1
) -> str:
    """
    Convert SH back to DWI signal.
    
    Parameters
    ----------
    sh_image : str
        Input SH image
    directions : str
        Gradient directions file
    output : str
        Output DWI image
    n_threads : int
        Number of threads
        
    Returns
    -------
    str
        Path to DWI image
    """
    run_mrtrix_cmd([
        "sh2amp", sh_image,
        directions,
        output,
        "-force",
        "-nthreads", str(n_threads)
    ])
    
    return output

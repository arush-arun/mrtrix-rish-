"""
QC Metrics Module

Compute quantitative metrics for harmonization quality assessment.
"""

from pathlib import Path
from typing import Dict, Optional
import subprocess


def compute_qc_metrics(
    original_sh: str,
    harmonized_sh: str,
    mask: Optional[str] = None,
    output_dir: Optional[str] = None
) -> Dict[str, float]:
    """
    Compute QC metrics comparing original and harmonized data.
    
    Parameters
    ----------
    original_sh : str
        Original SH image
    harmonized_sh : str
        Harmonized SH image
    mask : str, optional
        Brain mask
    output_dir : str, optional
        Directory for intermediate outputs
        
    Returns
    -------
    dict
        QC metrics
    """
    # TODO: Implement full metrics
    # For now, return placeholder
    
    metrics = {
        "status": "not_implemented",
        "fa_diff_mean": None,
        "fa_diff_std": None,
        "md_diff_mean": None,
        "md_diff_std": None,
        "angular_correlation": None,
    }
    
    return metrics


def compute_fa_md(sh_image: str, mask: str, output_dir: str) -> Dict[str, str]:
    """
    Compute FA and MD from SH coefficients.
    
    Parameters
    ----------
    sh_image : str
        SH coefficient image
    mask : str
        Brain mask
    output_dir : str
        Output directory
        
    Returns
    -------
    dict
        Paths to FA and MD images
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert SH to tensor
    tensor = output_dir / "tensor.mif"
    subprocess.run([
        "sh2tensor", sh_image, str(tensor),
        "-mask", mask, "-force"
    ], check=True)
    
    # Compute FA
    fa = output_dir / "fa.mif"
    subprocess.run([
        "tensor2metric", str(tensor),
        "-fa", str(fa), "-force"
    ], check=True)
    
    # Compute MD
    md = output_dir / "md.mif"
    subprocess.run([
        "tensor2metric", str(tensor),
        "-adc", str(md), "-force"
    ], check=True)
    
    return {
        "tensor": str(tensor),
        "fa": str(fa),
        "md": str(md)
    }


def compute_difference_maps(
    image1: str,
    image2: str,
    output: str,
    mask: Optional[str] = None
) -> str:
    """
    Compute voxel-wise difference between two images.
    
    Parameters
    ----------
    image1 : str
        First image
    image2 : str
        Second image
    output : str
        Output difference image
    mask : str, optional
        Brain mask
        
    Returns
    -------
    str
        Path to difference image
    """
    # Compute absolute difference
    subprocess.run([
        "mrcalc", image1, image2, "-sub", "-abs",
        output, "-force"
    ], check=True)
    
    # Apply mask if provided
    if mask:
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".mif", delete=False) as tmp:
            subprocess.run([
                "mrcalc", output, mask, "-mult",
                tmp.name, "-force"
            ], check=True)
            import shutil
            shutil.move(tmp.name, output)
    
    return output


def compute_stats(image: str, mask: Optional[str] = None) -> Dict[str, float]:
    """
    Compute statistics for an image.
    
    Parameters
    ----------
    image : str
        Input image
    mask : str, optional
        Brain mask
        
    Returns
    -------
    dict
        Statistics (mean, std, min, max)
    """
    cmd = ["mrstats", image]
    if mask:
        cmd.extend(["-mask", mask])
    
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    # Parse output
    # Format: "mean std min max count"
    lines = [l.strip() for l in result.stdout.strip().split("\n") if l.strip()]
    if len(lines) >= 2:
        values = lines[-1].split()
        return {
            "mean": float(values[0]),
            "std": float(values[1]),
            "min": float(values[2]),
            "max": float(values[3]),
            "count": int(float(values[4]))
        }
    
    return {}

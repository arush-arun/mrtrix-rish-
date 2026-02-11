"""
Harmonization Module

Applies RISH scale maps to spherical harmonic coefficients.
"""

import subprocess
from pathlib import Path
from typing import Dict, Optional, List
import tempfile
import shutil

from .rish_features import get_sh_indices, get_image_lmax, extract_rish_features
from .scale_maps import compute_scale_maps


def run_mrtrix_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run an MRtrix3 command."""
    return subprocess.run(cmd, capture_output=True, text=True, check=check)


def harmonize_sh(
    sh_image: str,
    scale_maps: Dict[int, str],
    output: str,
    lmax: Optional[int] = None,
    n_threads: int = 1
) -> str:
    """
    Apply RISH harmonization to SH coefficients.
    
    For each order l, multiply all m coefficients by scale_l.
    
    Parameters
    ----------
    sh_image : str
        Input SH coefficient image
    scale_maps : dict
        Mapping of order l -> scale map image path
    output : str
        Output harmonized SH image path
    lmax : int, optional
        Maximum SH order
    n_threads : int
        Number of threads
        
    Returns
    -------
    str
        Path to harmonized SH image
    """
    sh_image = Path(sh_image)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    
    if lmax is None:
        lmax = get_image_lmax(str(sh_image))
    
    sh_info = get_sh_indices(lmax)
    thread_opt = ["-nthreads", str(n_threads)] if n_threads > 1 else []
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        scaled_parts = []
        
        for l, (start_idx, end_idx) in sh_info.volume_indices.items():
            if l not in scale_maps:
                raise ValueError(f"No scale map for order {l}")
            
            scale_map = scale_maps[l]
            vol_range = f"{start_idx}:{end_idx}"
            
            # Extract coefficients for this order
            extracted = tmpdir / f"l{l}_coeffs.mif"
            run_mrtrix_cmd([
                "mrconvert", str(sh_image),
                "-coord", "3", vol_range,
                str(extracted),
                "-force"
            ] + thread_opt)
            
            # Apply scaling (broadcast scale map across m dimension)
            scaled = tmpdir / f"l{l}_scaled.mif"
            run_mrtrix_cmd([
                "mrcalc", str(extracted), scale_map, "-mult",
                str(scaled),
                "-force"
            ] + thread_opt)
            
            scaled_parts.append(str(scaled))
        
        # Concatenate all scaled parts
        run_mrtrix_cmd([
            "mrcat", *scaled_parts,
            "-axis", "3",
            str(output),
            "-force"
        ] + thread_opt)
    
    return str(output)


def harmonize_dwi(
    dwi_image: str,
    scale_maps: Dict[int, str],
    gradient_table: str,
    output: str,
    mask: Optional[str] = None,
    lmax: int = 8,
    n_threads: int = 1
) -> str:
    """
    Harmonize DWI data (convenience function).
    
    Converts DWI → SH, harmonizes, then optionally converts back.
    
    Parameters
    ----------
    dwi_image : str
        Input DWI image
    scale_maps : dict
        Scale maps from compute_scale_maps()
    gradient_table : str
        Path to gradient table (bvecs/bvals or MRtrix format)
    output : str
        Output path (SH or DWI depending on extension)
    mask : str, optional
        Brain mask
    lmax : int
        Maximum SH order for fitting
    n_threads : int
        Number of threads
        
    Returns
    -------
    str
        Path to harmonized output
    """
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    
    thread_opt = ["-nthreads", str(n_threads)] if n_threads > 1 else []
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Convert DWI to SH
        sh_image = tmpdir / "sh_coeffs.mif"
        cmd = [
            "amp2sh", dwi_image,
            "-lmax", str(lmax),
            str(sh_image),
            "-force"
        ]
        if mask:
            cmd.extend(["-mask", mask])
        run_mrtrix_cmd(cmd + thread_opt)
        
        # Harmonize SH
        sh_harmonized = tmpdir / "sh_harmonized.mif"
        harmonize_sh(
            str(sh_image),
            scale_maps,
            str(sh_harmonized),
            lmax=lmax,
            n_threads=n_threads
        )
        
        # If output should be SH, just copy
        if output.suffix in [".mif", ".nii", ".nii.gz"]:
            # Check if user wants SH or DWI output
            # For now, output SH
            shutil.copy(sh_harmonized, output)
        
    return str(output)


class RISHHarmonizer:
    """
    High-level interface for RISH harmonization.
    
    Examples
    --------
    >>> harmonizer = RISHHarmonizer(lmax=8)
    >>> harmonizer.create_template(
    ...     reference_dwi_list,
    ...     reference_mask_list,
    ...     "template/"
    ... )
    >>> harmonizer.harmonize(
    ...     "target/dwi.mif",
    ...     "target/mask.mif",
    ...     "harmonized/"
    ... )
    """
    
    def __init__(
        self,
        lmax: int = 8,
        smoothing_fwhm: float = 3.0,
        clip_range: tuple = (0.5, 2.0),
        n_threads: int = 4
    ):
        """
        Initialize RISH harmonizer.
        
        Parameters
        ----------
        lmax : int
            Maximum SH order
        smoothing_fwhm : float
            Scale map smoothing FWHM (mm)
        clip_range : tuple
            Scale factor clipping range
        n_threads : int
            Number of threads
        """
        self.lmax = lmax
        self.smoothing_fwhm = smoothing_fwhm
        self.clip_range = clip_range
        self.n_threads = n_threads
        
        self.template_dir = None
        self.scale_maps = None
        self.reference_rish = None
    
    def create_template(
        self,
        reference_sh_list: List[str],
        reference_mask_list: Optional[List[str]] = None,
        output_dir: str = "template"
    ) -> Dict[int, str]:
        """
        Create RISH template from reference site data.
        
        Parameters
        ----------
        reference_sh_list : list of str
            SH images from reference site (already in template space)
        reference_mask_list : list of str, optional
            Masks for each subject
        output_dir : str
            Output directory for template
            
        Returns
        -------
        dict
            Reference RISH features (averaged)
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        self.template_dir = output_dir
        
        # Extract RISH for each reference subject
        rish_list = []
        for i, sh_img in enumerate(reference_sh_list):
            mask = reference_mask_list[i] if reference_mask_list else None
            subj_dir = output_dir / "subjects" / f"sub-{i:03d}"
            rish = extract_rish_features(
                sh_img,
                str(subj_dir),
                lmax=self.lmax,
                mask=mask,
                n_threads=self.n_threads
            )
            rish_list.append(rish)
        
        # Average RISH features
        orders = list(rish_list[0].keys())
        self.reference_rish = {}
        
        for l in orders:
            images = [r[l] for r in rish_list]
            avg_output = output_dir / f"template_rish_l{l}.mif"
            
            run_mrtrix_cmd([
                "mrmath", *images,
                "mean", str(avg_output),
                "-force",
                "-nthreads", str(self.n_threads)
            ])
            
            self.reference_rish[l] = str(avg_output)
        
        return self.reference_rish
    
    def harmonize(
        self,
        target_sh: str,
        target_mask: Optional[str],
        output_dir: str
    ) -> Dict[str, str]:
        """
        Harmonize a target subject.
        
        Parameters
        ----------
        target_sh : str
            Target SH image (in template space)
        target_mask : str, optional
            Target brain mask
        output_dir : str
            Output directory
            
        Returns
        -------
        dict
            Paths to harmonized outputs
        """
        if self.reference_rish is None:
            raise RuntimeError("Must call create_template() first")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Extract target RISH
        target_rish_dir = output_dir / "rish"
        target_rish = extract_rish_features(
            target_sh,
            str(target_rish_dir),
            lmax=self.lmax,
            mask=target_mask,
            n_threads=self.n_threads
        )
        
        # Compute scale maps
        scale_maps_dir = output_dir / "scale_maps"
        scale_maps = compute_scale_maps(
            self.reference_rish,
            target_rish,
            str(scale_maps_dir),
            mask=target_mask,
            smoothing_fwhm=self.smoothing_fwhm,
            clip_range=self.clip_range,
            n_threads=self.n_threads
        )
        
        # Apply harmonization
        harmonized_sh = output_dir / "sh_harmonized.mif"
        harmonize_sh(
            target_sh,
            scale_maps,
            str(harmonized_sh),
            lmax=self.lmax,
            n_threads=self.n_threads
        )
        
        return {
            "harmonized_sh": str(harmonized_sh),
            "scale_maps": scale_maps,
            "target_rish": target_rish
        }
    
    def generate_qc_report(
        self,
        original: str,
        harmonized: str,
        output: str
    ) -> str:
        """
        Generate QC report comparing original and harmonized data.
        
        Parameters
        ----------
        original : str
            Original SH/DWI image
        harmonized : str
            Harmonized SH/DWI image
        output : str
            Output HTML report path
            
        Returns
        -------
        str
            Path to generated report
        """
        # TODO: Implement QC report generation
        raise NotImplementedError("QC report generation not yet implemented")


if __name__ == "__main__":
    print("Harmonization module. Import and use RISHHarmonizer class.")

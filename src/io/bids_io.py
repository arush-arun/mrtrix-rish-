"""
BIDS dataset I/O utilities.

Supports reading DWI data from BIDS-formatted datasets and
writing outputs to BIDS derivatives.
"""

from pathlib import Path
from typing import List, Dict, Optional, Tuple
import json
import subprocess


def find_bids_dwi(bids_dir: str, subject: Optional[str] = None) -> List[Dict]:
    """
    Find all DWI files in a BIDS dataset.
    
    Parameters
    ----------
    bids_dir : str
        Path to BIDS dataset root
    subject : str, optional
        Specific subject to find (e.g., "sub-01")
        
    Returns
    -------
    list of dict
        Each dict contains: subject, session, dwi, bval, bvec, json, mask
    """
    bids_dir = Path(bids_dir)
    results = []
    
    # Find subject directories
    if subject:
        subj_dirs = [bids_dir / subject] if (bids_dir / subject).exists() else []
    else:
        subj_dirs = sorted([d for d in bids_dir.iterdir() 
                          if d.is_dir() and d.name.startswith("sub-")])
    
    for subj_dir in subj_dirs:
        subject_id = subj_dir.name
        
        # Check for sessions
        session_dirs = sorted([d for d in subj_dir.iterdir() 
                              if d.is_dir() and d.name.startswith("ses-")])
        
        if not session_dirs:
            session_dirs = [subj_dir]  # No sessions, dwi directly under subject
        
        for ses_dir in session_dirs:
            session_id = ses_dir.name if ses_dir.name.startswith("ses-") else None
            
            # Find dwi directory
            if session_id:
                dwi_dir = ses_dir / "dwi"
            else:
                dwi_dir = subj_dir / "dwi"
            
            if not dwi_dir.exists():
                continue
            
            # Find DWI files (nii.gz or nii)
            # Prefer dir-AP over dir-PA (AP is typically main acquisition)
            for ext in [".nii.gz", ".nii"]:
                dwi_files = list(dwi_dir.glob(f"*_dwi{ext}"))
                
                # Filter: prefer dir-AP if both exist
                ap_files = [f for f in dwi_files if "dir-AP" in f.name]
                non_pe_files = [f for f in dwi_files if "dir-" not in f.name]
                
                # Use AP files if available, otherwise non-PE-encoded, otherwise all
                if ap_files:
                    dwi_files = ap_files
                elif non_pe_files:
                    dwi_files = non_pe_files
                # else keep all (including PA)
                
                for dwi_file in dwi_files:
                    entry = _parse_dwi_entry(dwi_file, subject_id, session_id)
                    if entry:
                        results.append(entry)
    
    return results


def _parse_dwi_entry(dwi_file: Path, subject: str, session: Optional[str]) -> Optional[Dict]:
    """Parse a single DWI file and find associated files."""
    stem = dwi_file.name.replace(".nii.gz", "").replace(".nii", "")
    dwi_dir = dwi_file.parent
    
    entry = {
        "subject": subject,
        "session": session,
        "dwi": str(dwi_file),
        "bval": None,
        "bvec": None,
        "json": None,
        "mask": None,
    }
    
    # Find bval
    bval_file = dwi_dir / f"{stem}.bval"
    if bval_file.exists():
        entry["bval"] = str(bval_file)
    
    # Find bvec
    bvec_file = dwi_dir / f"{stem}.bvec"
    if bvec_file.exists():
        entry["bvec"] = str(bvec_file)
    
    # Find JSON sidecar
    json_file = dwi_dir / f"{stem}.json"
    if json_file.exists():
        entry["json"] = str(json_file)
    
    # Check for mask in derivatives
    # Common locations: derivatives/masks/, derivatives/dwi_preproc/
    bids_root = dwi_dir.parent.parent
    if session:
        bids_root = bids_root.parent
    
    for deriv_name in ["masks", "dwi_preproc", "mrtrix", "qsiprep"]:
        deriv_dir = bids_root / "derivatives" / deriv_name
        if session:
            mask_dir = deriv_dir / subject / session / "dwi"
        else:
            mask_dir = deriv_dir / subject / "dwi"
        
        if mask_dir.exists():
            for mask_pattern in ["*_mask.nii*", "*_brain_mask.nii*", "*_brainmask.nii*"]:
                masks = list(mask_dir.glob(mask_pattern))
                if masks:
                    entry["mask"] = str(masks[0])
                    break
        if entry["mask"]:
            break
    
    return entry


def convert_to_mif(
    dwi_nifti: str,
    bval: str,
    bvec: str,
    output: str,
    json_file: Optional[str] = None
) -> str:
    """
    Convert BIDS NIfTI DWI to MRtrix .mif format.
    
    Parameters
    ----------
    dwi_nifti : str
        Input NIfTI file
    bval : str
        bval file
    bvec : str  
        bvec file
    output : str
        Output .mif file
    json_file : str, optional
        JSON sidecar for phase encoding info
        
    Returns
    -------
    str
        Path to output .mif file
    """
    cmd = [
        "mrconvert", dwi_nifti,
        "-fslgrad", bvec, bval,
        output, "-force"
    ]
    
    if json_file:
        cmd.extend(["-json_import", json_file])
    
    subprocess.run(cmd, check=True, capture_output=True)
    return output


def get_bids_derivative_path(
    bids_root: str,
    subject: str,
    filename: str,
    session: Optional[str] = None,
    pipeline: str = "mrtrix-rish",
    datatype: str = "dwi"
) -> Path:
    """
    Get the appropriate BIDS derivatives path for an output file.
    
    Parameters
    ----------
    bids_root : str
        BIDS dataset root
    subject : str
        Subject ID (e.g., "sub-01")
    filename : str
        Output filename
    session : str, optional
        Session ID (e.g., "ses-01")
    pipeline : str
        Pipeline name for derivatives folder
    datatype : str
        Data type (dwi, anat, etc.)
        
    Returns
    -------
    Path
        Full path for output file
    """
    deriv_root = Path(bids_root) / "derivatives" / pipeline
    
    if session:
        out_dir = deriv_root / subject / session / datatype
    else:
        out_dir = deriv_root / subject / datatype
    
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / filename


class BIDSDataset:
    """
    Simple BIDS dataset reader.
    
    Finds DWI data in BIDS-structured directories.
    """
    
    def __init__(self, root: str):
        """
        Initialize BIDS dataset.
        
        Parameters
        ----------
        root : str
            Path to BIDS dataset root
        """
        self.root = Path(root)
        
        if not self.root.exists():
            raise FileNotFoundError(f"Dataset not found: {root}")
        
        # Load dataset description if present
        desc_file = self.root / "dataset_description.json"
        self.description = {}
        if desc_file.exists():
            with open(desc_file) as f:
                self.description = json.load(f)
    
    def get_subjects(self) -> List[str]:
        """Get list of subject IDs."""
        subjects = []
        for path in self.root.iterdir():
            if path.is_dir() and path.name.startswith("sub-"):
                subjects.append(path.name)
        return sorted(subjects)
    
    def get_sessions(self, subject: str) -> List[str]:
        """Get list of sessions for a subject."""
        subj_dir = self.root / subject
        sessions = []
        
        for path in subj_dir.iterdir():
            if path.is_dir() and path.name.startswith("ses-"):
                sessions.append(path.name)
        
        return sorted(sessions) if sessions else [None]
    
    def get_dwi(
        self,
        subject: str,
        session: Optional[str] = None,
        derivatives: bool = False
    ) -> Optional[Dict[str, str]]:
        """
        Get DWI files for a subject.
        
        Parameters
        ----------
        subject : str
            Subject ID (e.g., "sub-01")
        session : str, optional
            Session ID (e.g., "ses-01")
        derivatives : bool
            Look in derivatives folder
            
        Returns
        -------
        dict or None
            Dictionary with keys: dwi, bval, bvec, json
        """
        if derivatives:
            base = self.root / "derivatives" / "mrtrix-rish"
        else:
            base = self.root
        
        # Build path
        if session:
            dwi_dir = base / subject / session / "dwi"
        else:
            dwi_dir = base / subject / "dwi"
        
        if not dwi_dir.exists():
            return None
        
        # Find DWI files
        result = {}
        
        # Look for nifti or mif
        for ext in [".nii.gz", ".nii", ".mif"]:
            dwi_files = list(dwi_dir.glob(f"*_dwi{ext}"))
            if dwi_files:
                result["dwi"] = str(dwi_files[0])
                break
        
        if "dwi" not in result:
            return None
        
        # Find bval/bvec
        stem = Path(result["dwi"]).stem.replace(".nii", "")
        
        for bval in dwi_dir.glob(f"{stem}*.bval"):
            result["bval"] = str(bval)
            break
        
        for bvec in dwi_dir.glob(f"{stem}*.bvec"):
            result["bvec"] = str(bvec)
            break
        
        # Find JSON sidecar
        for jsonf in dwi_dir.glob(f"{stem}*.json"):
            result["json"] = str(jsonf)
            break
        
        return result
    
    def get_all_dwi(
        self,
        derivatives: bool = False
    ) -> List[Dict]:
        """
        Get all DWI files in dataset.
        
        Returns
        -------
        list of dict
            Each dict has: subject, session, dwi, bval, bvec, json
        """
        results = []
        
        for subject in self.get_subjects():
            for session in self.get_sessions(subject):
                dwi = self.get_dwi(subject, session, derivatives)
                if dwi:
                    dwi["subject"] = subject
                    dwi["session"] = session
                    results.append(dwi)
        
        return results
    
    def write_derivative(
        self,
        content: bytes,
        subject: str,
        filename: str,
        session: Optional[str] = None,
        pipeline: str = "mrtrix-rish"
    ) -> str:
        """
        Write a file to derivatives.
        
        Parameters
        ----------
        content : bytes
            File content
        subject : str
            Subject ID
        filename : str
            Output filename
        session : str, optional
            Session ID
        pipeline : str
            Pipeline name for derivatives folder
            
        Returns
        -------
        str
            Path to written file
        """
        deriv_dir = self.root / "derivatives" / pipeline
        
        if session:
            out_dir = deriv_dir / subject / session / "dwi"
        else:
            out_dir = deriv_dir / subject / "dwi"
        
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / filename
        
        with open(out_path, "wb") as f:
            f.write(content)
        
        return str(out_path)

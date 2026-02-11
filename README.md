# mrtrix-rish

**MRtrix3-native RISH Harmonization Pipeline**

A reproducible, open-source tool for multi-site diffusion MRI harmonization using Rotational Invariant Spherical Harmonics (RISH), built entirely on the MRtrix3 ecosystem.

## Why?

Multi-site diffusion MRI studies suffer from scanner-induced variability. RISH harmonization operates at the spherical harmonics level, preserving angular (fiber orientation) information while removing site effects. Existing implementations (pnlbwh) rely on FSL/ANTs — this project provides a **pure MRtrix3 workflow**.

## Features

- 🧠 **FOD-level harmonization** — preserves fiber orientation information
- 🔧 **MRtrix3-native** — no FSL/ANTs dependency for core pipeline
- 📊 **Built-in QC** — automatic quality reports
- 📦 **BIDS-compatible** — works with standard neuroimaging data structures
- 🐳 **Containerized** — Docker/Singularity for reproducibility

## Quick Start

```bash
# Install
pip install mrtrix-rish

# Basic harmonization
mrtrix-rish harmonize \
    --reference site_a/ \
    --target site_b/ \
    --output harmonized/ \
    --template template/

# Or with config file
mrtrix-rish harmonize --config config/harmonize.yaml
```

## Requirements

- MRtrix3 >= 3.0.4
- Python >= 3.9
- NumPy, NiBabel (for QC only)

## Documentation

- [Design Document](docs/design/DESIGN.md)
- [Pipeline Flowchart](docs/design/FLOWCHART.md)
- [API Reference](docs/api/)
- [Examples](examples/)

## Citation

If you use this tool, please cite:
- Cetin Karayumak et al. (2019) - Original RISH method
- Tournier et al. (2019) - MRtrix3

## License

MIT License - see [LICENSE](LICENSE)

## Contributing

Contributions welcome! See [CONTRIBUTING.md](CONTRIBUTING.md)

# mrtrix-rish Docker image
# MRtrix3-native RISH harmonization for multi-site dMRI

FROM ubuntu:22.04

LABEL maintainer="Arush"
LABEL description="MRtrix3-native RISH harmonization pipeline"
LABEL version="0.1.0"

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    git \
    wget \
    curl \
    ca-certificates \
    libeigen3-dev \
    zlib1g-dev \
    libqt5opengl5-dev \
    libqt5svg5-dev \
    libgl1-mesa-dev \
    libfftw3-dev \
    libtiff5-dev \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Install MRtrix3 from conda-forge
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=/opt/micromamba/envs/mrtrix/bin:$PATH

RUN /bin/micromamba create -n mrtrix -c conda-forge mrtrix3 python=3.11 -y

# Set up Python environment
WORKDIR /app

# Copy package files
COPY pyproject.toml README.md ./
COPY src/ ./src/
COPY config/ ./config/

# Install mrtrix-rish
RUN /opt/micromamba/envs/mrtrix/bin/pip install -e ".[qc]"

# Create working directory
WORKDIR /data

# Default entrypoint
ENTRYPOINT ["mrtrix-rish"]
CMD ["--help"]

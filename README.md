1. Create & Activate Conda Environment (Recommended)
bash
conda create -n molecular_gpu python=3.9 -y
conda activate molecular_gpu
2. Install Core Dependencies (GPU-Accelerated)
bash
conda install -c conda-forge -y \
    numpy \
    pyqtgraph \
    pyqt \
    mdalibrary \
    cudatoolkit  # For NVIDIA GPUs (skip if using AMD)
3. Install Additional Libraries
bash
pip install \
    MDAnalysis \      # For .gro/.pdb file parsing
    PyOpenGL         # For OpenGL GPU rendering

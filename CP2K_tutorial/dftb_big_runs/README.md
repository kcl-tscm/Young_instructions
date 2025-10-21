# Introduction

This tutorial explains how to install CP2K using the provided set of commands, which are typical for managing software environments on a High-Performance Computing (HPC) cluster.

## CP2K Installation Tutorial (HPC Environment)

These commands are designed to prepare your environment, create a dedicated software space, and install the specific CP2K version you want.

| Command | Explanation |
| :--- | :--- |
| `module purge` | **Cleans Environment.** Removes any previously loaded software modules to prevent conflicts. |
| `module load python/miniconda3/24.3.0-0` | **Loads Conda.** Uses the HPC's module system to make the `conda` command available in your session. |
| `source $UCL_CONDA_PATH/etc/profile.d/conda.sh` | **Initializes Conda.** Activates necessary shell functions for `conda activate` to work. |
| `ls` | **List Files.** Standard command to list contents of the current directory (used as a check). |
| `conda create --name cp2k_env` | **Creates Environment.** Creates a new, isolated Conda environment named `cp2k_env`. |
| `conda activate cp2k_env` | **Activates Environment.** Switches the terminal session into the new `cp2k_env` environment. |
| `conda install -c conda-forge cp2k=2024.2=*openblas_openmpi*` | **Installs CP2K.** Installs the specified version of CP2K from the `conda-forge` channel, optimized with OpenBLAS and OpenMPI. |

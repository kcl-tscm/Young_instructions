# Introduction

This tutorial explains how to install CP2K using the provided set of commands, which are typical for managing software environments on a High-Performance Computing (HPC) cluster.

## CP2K Installation Tutorial (HPC Environment)

These commands are designed to prepare your environment, create a dedicated software space, and install the specific CP2K version you want.

## CP2K Installation Commands with Custom Path

This set of commands is tailored for HPC systems, using a module system and creating the Conda environment in a specific non-default directory (`/home/mmmXXXX/Scratch/conda_envs`).

| Command | Explanation |
| :--- | :--- |
| `module purge` | **Cleans Environment.** Removes any previously loaded software modules to avoid conflicts. |
| `module load python/miniconda3/24.3.0-0` | **Loads Conda.** Uses the module system to make the `conda` command available. |
| `source $UCL_CONDA_PATH/etc/profile.d/conda.sh` | **Initializes Conda.** Sets up the necessary shell functions for Conda environment management (like `conda activate`). |
| `ls` | **List Files.** A simple check to confirm the shell is responsive. |
| `conda create -p /home/mmmXXXX/Scratch/conda_envs/cp2k_env` | **Creates Environment.** Creates the isolated environment in your desired **Scratch** location using the `-p` (prefix) flag, naming the folder `cp2k_env`. |
| `conda activate /home/mmmXXXX/Scratch/conda_envs/cp2k_env` | **Activates Environment.** Activates the environment using its **full path** because it was created outside of the default Conda environment directory. |
| `conda install -c conda-forge cp2k=2024.2=*openblas_openmpi*` | **Installs CP2K.** Installs the specific version from the `conda-forge` channel, using the build optimized for **OpenBLAS** and **OpenMPI** for high performance. |

## How to Run CP2K

Once the installation is complete, follow these steps to use CP2K:

### 1. Activate the Environment

You must always activate the specific Conda environment where CP2K was installed before you can run the program.

```bash
conda activate cp2k_ant
```
### 2. Execute CP2K

Run the CP2K executable with your input file.

1. The -i flag specifies the input file.
2. The -o flag specifies the output file (where the results will be written).

```bash

cp2k -i your_input_file.inp -o your_output_file.out
```

Note: On an HPC system, you will typically submit this execution command within a batch script (PBS) to utilize multiple cores/nodes efficiently.
This is included in each of the created folders. 

### 3. Deactivate the Environment

When you are finished running your job or session, it is good practice to leave the dedicated environment.

```bash

conda deactivate
```
this is only needed for testing the executable. For the case of running in the scheduler, this command is not needed to be used.






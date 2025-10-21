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






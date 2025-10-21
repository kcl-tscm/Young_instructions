# ⚠️ Action Required: Update DFTB+ Input File (md.inp)

You are setting up a Molecular Dynamics (MD) simulation using DFTB+. To ensure the calculation runs successfully, you must update the file paths within your input file, md.inp.

## 1. The Requirement: Slater-Koster Files

The simulation needs access to the Slater-Koster (SK) files and dispersion parameters. These files are essential for the DFTB+ method.

## 2. File to Edit

You need to edit the file: **md.inp**

## 3. Required Changes

Locate the following section in your md.inp file and replace the existing paths with the new ones provided below.

The primary changes involve updating the absolute paths to your Young's username (mmmXXXX) and the final location of the SK files.

### Old Paths (Example):


```bash

DISPERSION_PARAMETER_FILE /home/mmm0666/Scratch/conda_envs/cp2k_env/share/cp2k/data/dftd3.dat
DISPERSION_TYPE D3
PARAM_FILE_PATH /home/mmm0666/Scratch/cp2k_anton_2024/3ob-3.1.0/skfiles
```

### New Paths (Mandatory Update):

```bash
DISPERSION_PARAMETER_FILE /home/mmmXXXX/Scratch/conda_envs/cp2k_env/share/cp2k/data/dftd3.dat
DISPERSION_TYPE D3
PARAM_FILE_PATH /home/mmmXXXX/Scratch/calculations/3ob-3.1.0/skfiles
```

## 🚨 Crucial Step

You must replace mmmXXXX in the new paths with your actual Young's username.

For example, if your username is mmm1234, the paths should look like this:

```bash

DISPERSION_PARAMETER_FILE /home/mmm1234/Scratch/conda_envs/cp2k_env/share/cp2k/data/dftd3.dat
DISPERSION_TYPE D3
PARAM_FILE_PATH /home/mmm1234/Scratch/calculations/3ob-3.1.0/skfiles
```

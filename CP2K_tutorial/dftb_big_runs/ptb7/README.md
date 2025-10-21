## WARNING

We are running a MD simulation using DFTB+. A set of files called slater-koster are needed for the calculation to be succesful. 
You need to change some variable paths in the md.inp file.

Instead of:

DISPERSION_PARAMETER_FILE /home/mmm0666/Scratch/conda_envs/cp2k_env/share/cp2k/data/dftd3.dat
DISPERSION_TYPE D3
PARAM_FILE_PATH /home/mmm0666/Scratch/cp2k_anton_2024/3ob-3.1.0/skfiles

You need to change it to:


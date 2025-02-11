**Installing LAMMPS with DeepMD and Other Libraries**

This guide will walk you through installing LAMMPS with the DeepMD library and several other essential packages. We'll also discuss how to adjust paths for personalized setups.

**Prerequisites**

* **Cluster:** This guide is designed for Young.
* **Module System:** Ensure you have a module system (e.g., `moduleload`) available.

**Instructions**

1. **Load Necessary Modules:**
   ```bash
   module purge
   module load beta-modules
   module load gcc-libs/10.2.0
   module load python/3.9.6-gnu-10.2.0
   module load emacs
   module load openblas/0.3.13-openmp/gnu-10.2.0
   module load compilers/gnu/10.2.0
   module load mpi/openmpi/4.0.5/gnu-10.2.0
   module load cmake/3.21.1
   module load fftw/3.3.9/
   ```

2. **Create and Activate Virtual Environment:**
   ```bash
   mkdir deepmdkit
   cd deepmdkit
   python3 -m venv envdeep
   source envdeep/bin/activate
   ```

3. **Adjust Path Variables:**
   ```bash
   # Replace '/home/mmm0666' with your home directory
   export PATH=$HOME/deepmdkit/envdeep/bin:$PATH
   export PYTHONPATH=$HOME/deepmdkit/envdeep/lib/python3.9/site-packages:$PYTHONPATH
   ```

4. **Install Libraries:**
   ```bash
   python3 -m pip install --upgrade pip
   export TF_ENABLE_ONEDNN_OPTS=0
   env MPICC=mpicc python3 -m pip install --force mpi4py
   python3 -m pip install tensorflow
   pip uninstall urllib3
   pip install "urllib3<2.0"
   pip install cython
   pip install numpy  
   pip install pybind11
   ```

5. **Clone DeepMD Repository:**
   ```bash
   git clone https://github.com/deepmodeling/deepmd-kit.git deepmd-kit
   ```

6. **Install DeepMD:**
   ```bash
   cd deepmd-kit
   python3 -m pip install .
   ```

7. **Build DeepMD Library for LAMMPS:**
   ```bash
   cd source/
   mkdir build
   cd build
   cmake -DENABLE_TENSORFLOW=TRUE -DUSE_TF_PYTHON_LIBS=TRUE -DCMAKE_INSTALL_PREFIX=/home/mmm0666/deepmdkit/envdeep ..
   make -j8
   make install
   make lammps
   cd ..
   cd ..
   cd ..
   ```
8. **Clone Plumed and install plumed as independent library**
   ```bash
      git clone https://github.com/plumed/plumed2.git
      cd plumed2
      ./configure --prefix=$HOME/new_deepmdkit/envdeep CXX=mpic++ CC=mpicc FC=mpif90 --enable-pycv
      make
      make install
   ```

10. **Clone LAMMPS Repository:**
   ```bash
   git clone -b release https://github.com/lammps/lammps.git mylammps
   ```

11. **Integrate DeepMD into LAMMPS:**
   ```bash
   cp -r deepmd-kit/source/build/USER-DEEPMD/ mylammps/src/
   cd mylammps/src/

   make yes-kspace
   make yes-extra-fix
   make yes-user-deepm
   make yes-plumed
   make lib-plumed args="-p $HOME/new_deepmdkit/envdeep"
   make mpi -j4
   ```

**Path Adjustments Explanation:**

*   In the instructions, we use `$HOME`. This ensures the paths are adjusted to your specific system.
*   This is crucial because the environment variables `PATH` and `PYTHONPATH` need to point to the correct locations of the installed libraries and executables.

**Important Note:**

*   If you encounter any errors during installation, double-check the module versions and path settings.
*   Refer to the documentation of DeepMD and LAMMPS for additional troubleshooting tips.

**Conclusion**

You've successfully installed LAMMPS with DeepMD and other essential libraries. This setup allows you to perform molecular dynamics simulations using the power of machine learning.

**Launching jobs in Young**

To launch jobs in Young, please use the submission file job.sh and change the paths accordingly where deep-MD kit and lmp_mpi are stored.




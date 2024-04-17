#!/bin/bash
#SBATCH --job-name=buildNWChemEx
#SBATCH --partition=kill-shared,shared
#SBATCH --time=0-18:00:00 ## time format is DD-HH:MM:SS
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=cn-02-03-01,cn-02-03-02,cn-02-03-03,cn-02-03-04,cn-02-03-05,cn-02-03-06,cn-02-03-07,cn-02-03-08,cn-02-03-09,cn-02-03-10,cn-02-08-01,cn-02-08-02,cn-02-08-03,cn-02-08-04,cn-02-08-05,cn-02-08-06,cn-02-08-07,cn-02-08-08,cn-02-08-09,cn-02-08-10,cn-02-13-01,cn-02-13-02,cn-02-13-03,cn-02-13-04,cn-02-13-05,cn-02-13-06,cn-02-13-07,cn-02-13-08,cn-02-13-09,cn-02-13-10,cn-03-03-01,cn-03-13-01,cn-03-13-02,cn-03-13-03,cn-03-13-04,cn-03-13-05,cn-03-13-06,cn-03-13-07,cn-03-13-08,cn-03-13-09,cn-03-13-10,cn-03-18-01,cn-03-18-02,cn-03-18-03,cn-03-18-04,cn-03-18-05,cn-03-18-06,cn-03-18-07,cn-03-18-08,cn-03-18-09,cn-03-18-10,cn-03-23-01,cn-03-23-02,cn-03-23-03,cn-03-23-04,cn-03-23-05,cn-03-23-06,cn-03-23-07,cn-03-23-08,cn-03-23-09,cn-03-23-10,cn-03-33-01,cn-03-33-02,cn-03-33-03,cn-03-33-04,cn-03-33-05,cn-03-33-06,cn-03-33-07,cn-03-33-08,cn-03-43-00,cn-04-21-00,cn-04-23-00,cn-04-35-00    # Only ivy-bridge CPUs
#SBATCH --nodes=1
#SBATCH --cores=8
#SBATCH --mem=320G
#SBATCH --distribution="*:*:*"  # Set the distribution to defaults if doing sbatch from interactive session
#SBATCH --error=build.err ## %A - filled with jobid
#SBATCH --output=build.out ## %A - filled with jobid



# We assume:
#   (1) we are buiding in the current directory
#   (2) the boost and libint directories are ready
buildDIR=$(pwd)

# Get the "toolchain.cmake" file ready (that was made in a previous build)
sed -i 's|NWX_MODULE_DIRECTORY.*|NWX_MODULE_DIRECTORY '"$buildDIR"'/NWChemEx-modules)|' ./toolchain.cmake


#####################################################################################

module purge
module load compiler/GCCcore/11.3.0 
module load mpi/OpenMPI/4.1.4-GCC-11.3.0
module load numlib/OpenBLAS/0.3.20-GCC-11.3.0
module load lang/Python/3.10.4-GCCcore-11.3.0 
module load devel/Boost/1.79.0-GCC-11.3.0
module load devel/CMake/3.24.3-GCCcore-11.3.0

### If this is the first time, create a virtual environment:
# python -m venv /mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/nwchemexTEST1/nwchemex_base/.nwchemexvenv
# source /mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/nwchemexTEST1/nwchemex_base/.nwchemexvenv/bin/activate
# pip install qcelemental
# pip install qcengine
# pip install networkx

### Otherwise, just load the virtual environment
source /mnt/lustre/koa/koastore/rsun_group/ruisun/sourcecode/.NWC_PY_V/bin/activate

#####################################################################################

# Clone the repo
git clone https://github.com/NWChemEx/NWChemEx.git 2>&1 | tee "OUTPUT.GITCLONE"
# <Type login information if prompted.>

# Navigate into the newly created NWChemEx subdirectory
cd NWChemEx

#     -DCMAKE_INSTALL_PREFIX=$buildDIR \
# Generate project buildsystem
cmake -H. \
      -Bbuild \
      -DCMAKE_TOOLCHAIN_FILE=`pwd`/../toolchain.cmake \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=$buildDIR/NWChemEx-install \
      2>&1 | tee "../OUTPUT.GEN"

# Kazuumi additions for the gradients:
cd build/_deps/
cd friendzone-src; git switch grad_mods; cd ..         # In friendzone
cd simde-src; git switch deriv_ptypes; cd ..           # In simde
cd ../..

# Build the project. You can change the "1" to another integer,
# N, to instead build with N threads
#cmake --build build -- -j 1 \
cmake --build build --parallel \
      2>&1 | tee "../OUTPUT.BUILD"

# Run tests
cd build && ctest 2>&1 | tee "../../OUTPUT.TEST"
#cd build && srun --overlap --mpi=pmix -n 4 ctest 2>&1 | tee "../../OUTPUT.TEST"

# Back out of the build directory
cd ..

# Install the project
cmake --build build --target install 2>&1 | tee "../OUTPUT.INSTALL"

# Return to the top level directory
cd ..

#!/bin/bash
#SBATCH --job-name=NWChemExMD
#SBATCH --partition=shared,kill-shared
##SBATCH --exclusive
##SBATCH --nodelist=cn-03-33-09,cn-03-33-10,cn-03-38-01,cn-03-38-02,cn-03-38-03,cn-03-38-04,cn-03-38-05,cn-03-38-06,cn-03-38-07,cn-03-38-08,cn-03-38-09,cn-03-38-10,cn-05-03-01,cn-05-03-02,cn-05-03-03,cn-05-03-04,cn-05-03-05,cn-05-03-06,cn-05-03-07,cn-05-03-08,cn-05-03-09,cn-05-03-10,cn-05-08-01,cn-05-08-02,cn-05-08-03,cn-05-08-04,cn-05-08-05,cn-05-08-06,cn-05-08-07,cn-05-08-08,cn-05-08-09,cn-05-08-10,cn-05-13-01,cn-05-13-02,cn-05-13-03,cn-05-13-04,cn-05-13-05,cn-05-13-06,cn-05-13-07,cn-05-13-08,cn-05-13-09,cn-05-13-10,cn-06-03-01,cn-06-03-02,cn-06-03-03,cn-06-03-04,cn-06-03-05,cn-06-03-06,cn-06-03-07,cn-06-03-08,cn-06-03-09,cn-06-03-10,cn-06-08-01,cn-06-08-02,cn-06-08-03,cn-06-08-04,cn-06-08-05,cn-06-08-06,cn-06-08-07,cn-06-08-08,cn-06-08-09,cn-06-08-10,cn-06-13-01,cn-06-13-02,cn-06-13-03,cn-06-13-04,cn-06-13-05,cn-06-13-06,cn-06-13-07,cn-06-13-08,cn-06-13-09,cn-06-13-10,cn-07-03-01,cn-07-03-02,cn-07-03-03,cn-07-03-04,cn-07-03-05,cn-07-03-06,cn-07-03-07,cn-07-03-08,cn-07-03-09,cn-07-03-10,cn-07-18-01,cn-07-18-02,cn-07-18-03,cn-07-27-01,cn-07-28-01,cn-07-28-02,cn-07-29-01,cn-07-29-02,cn-07-30-01,cn-07-30-02,cn-07-31-01,cn-08-30-01,cn-08-30-03,cn-08-30-04,cn-08-30-05,cn-08-30-06,cn-08-30-07,cn-08-30-08,cn-08-30-09,cn-08-30-10   # Only haswell architectures
#SBATCH --nodelist=cn-03-33-09,cn-03-33-10,cn-03-38-01,cn-03-38-02,cn-03-38-03,cn-03-38-04,cn-03-38-05,cn-03-38-06,cn-03-38-07,cn-03-38-08,cn-03-38-09,cn-03-38-10,cn-05-03-01,cn-05-03-02,cn-05-03-03,cn-05-03-04,cn-05-03-05,cn-05-03-06,cn-05-03-07,cn-05-03-08,cn-05-03-09,cn-05-03-10,cn-05-08-01,cn-05-08-02,cn-05-08-03,cn-05-08-04,cn-05-08-05,cn-05-08-06,cn-05-08-07,cn-05-08-08,cn-05-08-09,cn-05-08-10,cn-05-13-01,cn-05-13-02,cn-05-13-03,cn-05-13-04,cn-05-13-05,cn-05-13-06,cn-05-13-07,cn-05-13-08,cn-05-13-09,cn-05-13-10,cn-06-03-01,cn-06-03-02,cn-06-03-03,cn-06-03-04,cn-06-03-05,cn-06-03-06,cn-06-03-07,cn-06-03-08,cn-06-03-09,cn-06-03-10,cn-06-08-01,cn-06-08-02,cn-06-08-03,cn-06-08-04,cn-06-08-05,cn-06-08-06,cn-06-08-07,cn-06-08-08,cn-06-08-09,cn-06-08-10,cn-06-13-01,cn-06-13-02,cn-06-13-03,cn-06-13-04,cn-06-13-05,cn-06-13-06,cn-06-13-07,cn-06-13-08,cn-06-13-09,cn-06-13-10,cn-07-03-01,cn-07-03-02,cn-07-03-03,cn-07-03-04,cn-07-03-05,cn-07-03-06,cn-07-03-07,cn-07-03-08,cn-07-03-09,cn-07-03-10,cn-07-18-01,cn-07-18-02,cn-07-18-03,cn-07-27-01,cn-07-28-01,cn-07-28-02,cn-07-29-01,cn-07-29-02,cn-07-30-01,cn-07-30-02,cn-07-31-01,cn-08-30-01,cn-08-30-03,cn-08-30-04,cn-08-30-05,cn-08-30-06,cn-08-30-07,cn-08-30-08,cn-08-30-09,cn-08-30-10,cn-09-03-01,cn-09-03-02,cn-09-04-01,cn-09-04-02,cn-09-05-01,cn-09-05-02,cn-09-06-01,cn-09-06-02,cn-09-07-02,cn-09-30-01,cn-09-30-02,cn-09-31-01,cn-09-31-02,cn-09-32-01,cn-09-32-02,cn-09-33-01,cn-09-33-02,cn-15-03-01,cn-15-03-02,cn-15-04-01,cn-15-04-02,cn-15-05-01,cn-15-05-02,cn-15-06-01,cn-15-06-02,cn-15-07-01,cn-15-07-02,cn-15-08-01,cn-15-08-02,cn-15-09-01,cn-15-09-02,cn-15-10-01,cn-15-10-02,cn-15-11-01,cn-15-11-02,cn-15-12-01,cn-15-12-02,cn-15-13-01,cn-15-13-02,cn-15-14-01,cn-15-14-02,cn-15-15-01,cn-15-15-02,cn-15-16-01,cn-15-16-02,cn-15-21-01,cn-15-21-02,cn-15-22-01,cn-15-22-02,cn-15-23-01,cn-15-23-02,cn-15-24-01,cn-15-24-02,cn-15-25-01,cn-15-25-02,cn-15-26-01,cn-15-26-02,cn-15-27-01,cn-15-27-02,cn-15-28-01,cn-15-28-02,cn-15-29-02,cn-15-31-01,cn-15-31-02,cn-15-32-01,cn-15-32-02,cn-15-33-02,cn-15-34-02   # Only non-ivy-bridge CPU architectures
##SBATCH --partition=rsun
##SBATCH --account=rsun
#SBATCH --no-requeue
#SBATCH --mem=125000
#SBATCH --nodes=1
#SBATCH --cores=9
##SBATCH --tasks-per-node=1
##SBATCH --cpus-per-task=1
#SBATCH --error=err
##SBATCH --exclude=cn-09-04-01
#SBATCH --time=1-23:59:59 ## time format is DD-HH:MM:SS
#SBATCH --output=out

if [ "$#" -ne "1" ]; then
  echo "ERROR: Wrong number of arguments. Must be:"
  echo "   (1) Impact parameter in angstrom"
  exit
fi




unset SLURM_CPU_BIND_VERBOSE SLURM_CPU_BIND_LIST SLURM_CPU_BIND_TYPE SLURM_CPU_BIND
ulimit -s unlimited
source ~/.bash_profile

module purge
module load compiler/GCCcore/11.3.0 
module load mpi/OpenMPI/4.1.4-GCC-11.3.0
module load numlib/OpenBLAS/0.3.20-GCC-11.3.0
module load lang/Python/3.10.4-GCCcore-11.3.0 
module load devel/Boost/1.79.0-GCC-11.3.0
module load devel/CMake/3.24.3-GCCcore-11.3.0

source /mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/MLexp1/MLMD-cli/.nwchemexASEvenv/bin/activate

#moduleDIR=/mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/nwchemexTEST1/nwchemex_gcc_gccNWC_automated_try1
moduleDIR=/mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/nwchemexTEST1/nwchemex_gcc_gccNWC_automated_try7

export PYTHONPATH="$PYTHONPATH:$moduleDIR/NWChemEx-modules"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$moduleDIR/NWChemEx-install/lib64"
while read dir; do
  echo "Adding $dir to LD_LIBRARY_PATH..."
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$dir"
done < <(ls -d $moduleDIR/NWChemEx-install/lib64/*/ | sed 's|/$||')

# NWChem paths:
export NWCHEM_TOP="/mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/nwchemGCCtest1/nwchem-ivybridge-try1/nwchem"
export NWCHEM_EXECUTABLE=$NWCHEM_TOP/bin/LINUX64/nwchem
export NWCHEM_BASIS_LIBRARY=$NWCHEM_TOP/src/basis/libraries/

export PATH="$PATH:$NWCHEM_TOP/bin/LINUX64"

# SLURM variables:
export PMIX_MCA_psec=^munge

module list

nodename=$(hostname)
echo "Node: $nodename"
grep "$nodename" /mnt/lustre/koa/koastore/rsun_group/nodeActiveFeaturesList.txt

# Load your psi4 environment (if using it)
# source activate psi4

b=$1

# Assume that we are only doing one task per job on a single node
echo "Number of tasks: ${SLURM_NTASKS}"
let "Ntasks = ${SLURM_NTASKS} - 1"

#############################################################

xyzfile=/mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/MLexp1/MLMD-cli/tmp.A.xyz
PESfile=/mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/MLexp1/MLMD-cli/tmp.A.nwchemex

echo "START: "$(date)

rCM=5.0
CE=4.9713
#b=1.0

dt=0.15
maxSteps=100000

# Start the MD by doing some initial sampling
#srun --overlap --mpi=pmix -n $Ntasks python -u /home/kazuumi/rsun_koastore/kazuumiTest1/MLexp1/MLMD-cli/cli.py $xyzfile $PESfile . --atomsInFirstGroup "1 2" --collisionEnergy $CE --impactParameter $b --centerOfMassDistance $rCM --production $maxSteps --interval 1 --time_step $dt --INITQPa "semiclassical" --INITQPb "thermal" --NVIBa 0 --NROTa 0 --TVIBb 10.0 --TROTb 10.0 --n_threads $Ntasks > production.log
srun --overlap --mpi=pmix -n $Ntasks python -u /home/kazuumi/rsun_koastore/kazuumiTest1/MLexp1/MLMD-cli/cli.py $xyzfile $PESfile . --atomsInFirstGroup "1 2" --collisionEnergy $CE --impactParameter $b --centerOfMassDistance $rCM --production $maxSteps --interval 1 --time_step $dt --INITQPa "thermal" --INITQPb "thermal" --TVIBa 298.15 --TROTa 298.15 --TVIBb 298.15 --TROTb 298.15 --n_threads $Ntasks > production.log

echo "FINISH: "$(date)

python -u /mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/pythonSELECTandMD/psi4traj2exyz.py ./production.traj > ./production.xyz

echo "XYZ: "$(date)


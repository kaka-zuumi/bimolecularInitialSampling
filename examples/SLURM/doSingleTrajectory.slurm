#!/bin/bash
#SBATCH --job-name=MLMDsGDML
#SBATCH --partition=shared,kill-shared
##SBATCH --partition=rsun
##SBATCH --account=rsun
#SBATCH --no-requeue
#SBATCH --mem=25000
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --error=err
#SBATCH --exclude=cn-09-04-01
#SBATCH --time=1-23:59:59 ## time format is DD-HH:MM:SS
#SBATCH --output=out

if [ "$#" -ne "1" ]; then
  echo "ERROR: Wrong number of arguments. Must be:"
  echo "   (1) Impact parameter in angstrom"
  exit
fi

module purge
module load lang/Anaconda3

# Load your psi4 environment (if using it)
# source activate psi4

b=$1

# Assume that we are only doing one task per job on a single node
n_threads=$SLURM_CPUS_PER_TASK

#############################################################

# Site with all methods and basis sets:
# https://psicode.org/psi4manual/master/dft_byfunctional.html#table-dft-meta
# https://psicode.org/psi4manual/master/basissets_tables.html

xyzfile=/home/kazuumi/rsun_koastore/kazuumiTest1/MLexp1/MLMD-cli/input.xyz
PESfile=/home/kazuumi/rsun_koastore/kazuumiTest1/MLexp1/globalPES/PES-ntrain14000-try1/sig15.npz

echo "START: "$(date)

rCM=10.0
CE=4.349904
#b=1.0

dt=0.25
maxSteps=100000

# Start the MD by doing some initial sampling
python -u /home/kazuumi/rsun_koastore/kazuumiTest1/MLexp1/MLMD-cli/cli.py $xyzfile $PESfile . --atomsInFirstGroup "1 2" --collisionEnergy $CE --impactParameter $b --centerOfMassDistance $rCM --production $maxSteps --interval 1 --time_step $dt --INITQPa "semiclassical" --INITQPb "thermal" --NVIBa 0 --NROTa 0 --TVIBb 10.0 --TROTb 10.0 --n_threads $n_threads > production.log

echo "FINISH: "$(date)

python -u /home/kazuumi/rsun_koastore/kazuumiTest1/MLexp1/traj2exyz.py ./production.traj > ./production.xyz

echo "XYZ: "$(date)

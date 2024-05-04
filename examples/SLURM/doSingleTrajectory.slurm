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

# Load your psi4 environment (or other package)
# source activate psi4
# source activate ../../.bisampleschnet
# source activate ../../.bisamplesgdml

b=$1

# Assume that we are only doing one task per job on a single node
n_threads=$SLURM_CPUS_PER_TASK

#############################################################

# Site with all methods and basis sets:
# https://psicode.org/psi4manual/master/dft_byfunctional.html#table-dft-meta
# https://psicode.org/psi4manual/master/basissets_tables.html

xyzfile=../../examples/input.hbr.ch4.xyz
#PESfile=../../examples/input.BrCH5.psi4
PESfile=../../examples/schnet.BrCH5.best_model

#xyzfile=../../examples/input.ch.h2s.xyz
#PESfile=../../examples/sgdml.ch.h2s.npz
#PESfile=../../examples/schnet.ch.h2s.best_model

echo "START: "$(date)

rCM=10.0
CE=4.349904
#b=1.0

dt=0.25
maxSteps=100000

# Start the MD by doing some initial sampling
python -u ../../cli.py $xyzfile $PESfile . --atomsInFirstGroup "1 2" --collisionEnergy $CE --impactParameter $b --centerOfMassDistance $rCM --production $maxSteps --interval 1 --time_step $dt --INITQPa "semiclassical" --INITQPb "thermal" --NVIBa 0 --NROTa 0 --TVIBb 10.0 --TROTb 10.0 --n_threads $n_threads > production.log

echo "FINISH: "$(date)


#!/bin/bash
#SBATCH --job-name=NWChemExMD
#SBATCH --partition=shared,kill-shared
##SBATCH --partition=rsun
##SBATCH --account=rsun
#SBATCH --no-requeue
#SBATCH --mem=125000
#SBATCH --nodes=1
#SBATCH --cores=8     # For GAMESS with mpirun, must use an even number
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
module load lang/Anaconda3
module load toolchain/intel/2022b
module load toolchain/iimpi/2022b

module load lang/Python/3.10.4-GCCcore-11.3.0

source /mnt/lustre/koa/koastore/rsun_group/kazuumiTest1/MLexp1/MLMD-cli/.nwchemexASEvenv/bin/activate

# GAMESS paths:
export GAMESS_TOP="/mnt/lustre/koa/koastore/rsun_group/camels/gamess"
export PATH="$PATH:$GAMESS_TOP"

module list

#####################################################################################

nodename=$(hostname)
echo "Node: $nodename"
grep "$nodename" /mnt/lustre/koa/koastore/rsun_group/nodeActiveFeaturesList.txt

b=$1

# Assume that we are only doing one task per job on a single node
echo "Number of tasks: ${SLURM_NTASKS}"
Ntasks=${SLURM_NTASKS}

#############################################################

# The PES file specifies (1) the reference wavefunction method (referencemethod=rhf,uhf,rohf),
# (2) the SCF method (e.g., method=SCF, B3LYP, MP2, ...), and (3) the basis set (e.g. basis_set=MINI,DZV,...)
# GAMESS basis set documentation: https://myweb.liu.edu/~nmatsuna/gamess/input/BASIS.html
# QCEngine/GAMESS methods documentation: https://github.com/MolSSI/QCEngine/blob/f5f6da3751373fa9b57ea484cbf71416ba679743/qcengine/programs/gamess/germinate.py
# Note: QCEngine/GAMESS can't specify the gaussian basis sets (6-31G) easily without future modification
xyzfile=../../examples/input.ho.ch2.xyz
PESfile=../../examples/input.CH3O.gamess.qcengine

echo "START: "$(date)

rCM=5.0
CE=4.9713
#b=1.0

dt=0.15
maxSteps=100000

# Start the MD by doing some initial sampling
#python -u ../../cli.py $xyzfile $PESfile . --atomsInFirstGroup "1 2" --collisionEnergy $CE --impactParameter $b --centerOfMassDistance $rCM --production $maxSteps --interval 1 --time_step $dt --INITQPa "semiclassical" --INITQPb "thermal" --NVIBa 0 --NROTa 0 --TVIBb 10.0 --TROTb 10.0 --n_threads $Ntasks > production.log
python -u ../../cli.py $xyzfile $PESfile . --atomsInFirstGroup "1 2" --collisionEnergy $CE --impactParameter $b --centerOfMassDistance $rCM --production $maxSteps --interval 1 --time_step $dt --INITQPa "thermal" --INITQPb "thermal" --TVIBa 298.15 --TROTa 298.15 --TVIBb 298.15 --TROTb 298.15 --n_threads $Ntasks > production.log

echo "FINISH: "$(date)


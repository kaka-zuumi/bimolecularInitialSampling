
## Installing the package

Download the git repository and enter the folder.
```                                                                                                                                                                      
git clone https://github.com/kaka-zuumi/bimolecularInitialSampling.git                                                               
cd bimolecularInitialSampling
```

The initial sampling package makes use of energy and force calls from an ASE `atoms` object with an attached `calculator`. The same calculator can then be used to do molecular dynamics with. The calculator can either be an analytical functional (e.g. some ML potential) or an _ab initio_ electronic structure calculation. As of now, psi4, NWChemEx, and GAMESS (via QCEngine) interfaces are available for use. The examples here make use of the `psi4` python package for the latter case.

For the psi4 interface, we suggest creating a conda environment with two packages and activating it like so:
```
conda create -n bisamplepsi4ase psi4 ase                                                               
conda activate bisamplepsi4ase                                                                                   
```

For the Schnet interface, we suggest creating a python virtual environment with the schnetpack package and activating it like so:
```
python -m venv .bisampleschnet
source .bisampleschnet/bin/activate
pip install schnetpack==1.0
```


For the QCEngine/GAMESS interface (which requires GAMESS pre-installed), we suggest creating a python virtual environment with four packages and activating it like so:
```
python -m venv .bisampleqcenginegamess
source .bisampleqcenginegamess/bin/activate
pip install qcengine qcelemental networkx ase
```

For the NWChemEx interface, there is quite a lot of work involved in building the package... After finding the location of appropriate packages and modules, the script "examples/SLURM/build_nwx.slurm" and toolchain file "examples/SLURM/toolchain.cmake" can be changed, and then (with both in a fresh directory) the build_nwx.slurm script can be executed/submitted to do the build. Depending on the number of processes, it may take anywhere from 1-6 hours. This assumes NWChem is already installed. A more detailed install guide is in: https://nwchemex.github.io/NWChemEx/installation/building.html


## Try it out yourself!

To do an example initial sampling with molecular dynamics using a psi4 potential energy surface, try one of the exampler reactions below:

###  Br + CH<sub>5</sub><sup>+</sup>  ⟶

```
mkdir test1/; cd test1/
python -u ../cli.py ../examples/input.br.ch5.xyz ../examples/input.BrCH5.psi4 . --atomsInFirstGroup "1" --collisionEnergy 5.0 --impactParameter 1.0 --centerOfMassDistance 10.0 --production 1000 --interval 1 --time_step 0.15 --INITQPa "thermal" --INITQPb "thermal" --TVIBa 298.15 --TROTa 298.15 --TVIBb 298.15 --TROTb 298.15 > asepsi4md0.out
```

###  HBr<sup>+</sup> + CH<sub>4</sub>  ⟶

```
mkdir test2/; cd test2/
python -u ../cli.py ../examples/input.hbr.ch4.xyz ../examples/input.BrCH5.psi4 . --atomsInFirstGroup "1 2" --collisionEnergy 5.0 --impactParameter 1.0 --centerOfMassDistance 10.0 --production 1000 --interval 1 --time_step 0.15 --INITQPa "semiclassical" --INITQPb "thermal" --NVIBa 0 --NROTa 3 --TVIBb 298.15 --TROTb 298.15 > asepsi4md0.out
```

###  H<sub>2</sub>Br<sup>+</sup> + CH<sub>3</sub>  ⟶

```
mkdir test3/; cd test3/
python -u ../cli.py ../examples/input.h2br.ch3.xyz ../examples/input.BrCH5.psi4 . --atomsInFirstGroup "1 2 4" --collisionEnergy 5.0 --impactParameter 1.0 --centerOfMassDistance 10.0 --production 1000 --interval 1 --time_step 0.15 --INITQPa "thermal" --INITQPb "thermal" --TVIBa 298.15 --TROTa 298.15 --TVIBb 298.15 --TROTb 298.15 > asepsi4md0.out
```


Command line arguments can be explained with "cli.py --help". In general, three positional arguments are always required: (1) the combined XYZ file of both reactants, (2) the potential energy surface file with a specific file ending depending on the method (".nwchemex",".gamess.qcengine",".psi4", and ".npz" correspond to the NWChemEx, QCEngine/GAMESS, psi4, and sGDML interfaces, respectively), and (3) the directory to place output files like the trajectory. 




## Simulations on an HPC cluster

While on an interactive terminal, these simulations can be done one-by-one, when submitted to a node on a cluster, hundreds of thousands of these simulations can be done at once. After getting thousands of trajectories, statistically meaningful averages (e.g., product yields, intermediate lifetimes, rate constants) can be calculated. However, configuring things on your own HPC cluster (with its own job scheduler) may be tricky. Examples shown below are for an HPC cluster with a SLURM job scheduler.

Simulations making use of an sGDML or psi4 potential energy surface only require loading the python package, either with the conda or virtual environment described earlier. An example trajectory can be submitted with "examples/SLURM/doSingleTrajectory.slurm" using a varying impact paramter, here 1.0, with:
```
sbatch examples/SLURM/doSingleTrajectory.slurm 1.0
```

Simulations making use of a QCEngine/GAMESS potential energy surface require both loading the python package, as well as specifying PATHs in the environment (with an "export" statement in bash) for the GAMESS executable, as well as changing PATHs in the "rungms.MPI" executable so as to use an appropriate scratch directory (in this case, the current directory "."). An exaple trajectory can be submitted with "examples/SLURM/doQCEngineGAMESSSingleTrajectory.slurm" using a varying impact parameter, here 1.0, with:
```
sbatch examples/SLURM/doQCEngineGAMESSSingleTrajectory.slurm 1.0
```

Simulations making use of a NWChemEx potential energy surface require the loading the python package used in the build script "build_nwx.slurm" as well as loading all relevant modules and PATHs as used in the build script. Additional paths must be specified for the new install and module directories made during the build. An example trajectory can be submitted with "examples/SLURM/doNWChemExSingleTrajectory.slurm" using a varying impact parameter, here 1.0, with:
```
sbatch examples/SLURM/doNWChemExSingleTrajectory.slurm 1.0
```



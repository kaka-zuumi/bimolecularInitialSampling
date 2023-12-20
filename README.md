
## Installing the package

The initial sampling package makes use of energy and force calls from an ASE `atoms` object with an attached `calculator`. The same calculator can then be used to do molecular dynamics with. The calculator can either be an analytical functional (e.g. some ML potential) or an _ab initio_ electronic structure calculation. The examples here make use of the `psi4` python package for the latter case.

For the ASE/psi4 interface, we suggest creating a conda environment with the two packages like so:



## Try it out yourself!

To do an example initial sampling with molecular dynamics using a psi4 potential energy surface, try one of the exampler reactions below:

###  Br + CH<sub>5</sub><sup>+</sup>  ⟶

```
mkdir test1/; cd test1/
python -u ../cli.py ../examples/input.br.ch5.xyz ../examples/input.psi4 . --atomsInFirstGroup "1" --collisionEnergy 5.0 --impactParameter 1.0 --centerOfMassDistance 10.0 --production 1000 --interval 1 --time_step 0.15 --INITQPa "thermal" --INITQPb "thermal" --TVIBa 298.15 --TROTa 298.15 --TVIBb 298.15 --TROTb 298.15 > asepsi4md0.out
```

###  HBr<sup>+</sup> + CH<sub>4</sub>  ⟶

```
mkdir test2/; cd test2/
python -u ../cli.py ../examples/input.hbr.ch4.xyz ../examples/input.psi4 . --atomsInFirstGroup "1 2" --collisionEnergy 5.0 --impactParameter 1.0 --centerOfMassDistance 10.0 --production 1000 --interval 1 --time_step 0.15 --INITQPa "semiclassical" --INITQPb "thermal" --NVIBa 0 --NROTa 3 --TVIBb 298.15 --TROTb 298.15 > asepsi4md0.out
```

###  H<sub>2</sub>Br<sup>+</sup> + CH<sub>3</sub>  ⟶

```
mkdir test3/; cd test3/
python -u ../cli.py ../examples/input.h2br.ch3.xyz ../examples/input.psi4 . --atomsInFirstGroup "1 2 3" --collisionEnergy 5.0 --impactParameter 1.0 --centerOfMassDistance 10.0 --production 1000 --interval 1 --time_step 0.15 --INITQPa "thermal" --INITQPb "thermal" --TVIBa 298.15 --TROTa 298.15 --TVIBb 298.15 --TROTb 298.15 > asepsi4md0.out
```


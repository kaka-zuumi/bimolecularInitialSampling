

To do an example initial sampling with molecular dynamics using a psi4 potential energy surface, try one of the exampler reactions below:

#  HBr^+ + CH_4

```
python -u cli.py input.xyz input.psi4 . --atomsInFirstGroup "1 2" --collisionEnergy 5.0 --impactParameter 1.0 --centerOfMassDistance 10.0 --production 1000 --interval 1 --time_step 0.15 --INITQPa "QM" --INITQPb "thermal" --NVIBa 0 --NROTa 3 --TVIBb 298.15 --TROTb 298.15 > asepsi4md0.out
```

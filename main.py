from ase.io import read
from ase.md.verlet import VelocityVerlet
from ase import units

from ase.vibrations import Vibrations
import random

import argparse
import numpy as np
import os

from psi4calc import psi4calculator
from initialSampling import initialSampling

###################################################

# Define global constants up here in the correct
# units so that internally, everything uses:
#    Energy: eV
#  Distance: Angstrom
#      Mass: Dalton

global r2threshold

###################################################

parser = argparse.ArgumentParser(description="Do a single MD trajectory using a initial geometry (and momenta) and a sGDML model",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("initialGeometryFile", type=str, help="XYZ file with initial geometry; if initial conditions are sampled in the script, then this argument is required but is just an example XYZ")
parser.add_argument("psi4inputFile", type=str, help="psi4 input file")
parser.add_argument("outputDir", type=str, help="Directory to output stuff in")
parser.add_argument("--isotopeMassesFile", type=str, help="Change masses of specific atoms e.g. like isotopic substitution", default=None)
parser.add_argument("--initialMomentaFile", type=str, help="XYZ file with initial momenta")
parser.add_argument("--atomsInFirstGroup", type=str, help="String with atoms which are in first group of atoms, separated by spaces")
parser.add_argument("--collisionEnergy", type=float, help="Collision energy in kcal/mol")
parser.add_argument("--impactParameter", type=float, help="Impact parameter in Angstrom")
parser.add_argument("--centerOfMassDistance", type=float, help="Distance between the two molecules' centers of mass in Angstrom")
parser.add_argument("--production", type=int, help="Supply number of steps for a production run")
parser.add_argument("--interval", type=int, help="How often to print out the energy")
parser.add_argument("--time_step", type=float, help="The time step in fs for a production run")
parser.add_argument("--n_threads", type=int, help="The number of threads to ask psi4 to use")

parser.add_argument("--INITQPa", type=str, help="Initial sampling method for atoms in first group ('semiclassical', 'thermal', or None)", default=None)
parser.add_argument("--NVIBa", type=int, help="Vibrational quantum number of atoms in first group (supply if using the 'QM' initial sampling)")
parser.add_argument("--NROTa", type=int, help="Rotational quantum number of atoms in first group (supply if using the 'QM' initial sampling)")
parser.add_argument("--TVIBa", type=float, help="Vibrational temperature of atoms in first group (supply if using the 'thermal' initial sampling)")
parser.add_argument("--TROTa", type=float, help="Rotational temperature of atoms in first group (supply if using the 'thermal' initial sampling)")

parser.add_argument("--INITQPb", type=str, help="Initial sampling method for atoms in second group ('semiclassical', 'thermal', or None)", default=None)
parser.add_argument("--NVIBb", type=int, help="Vibrational quantum number of atoms in second group (supply if using the 'QM' initial sampling)")
parser.add_argument("--NROTb", type=int, help="Rotational quantum number of atoms in second group (supply if using the 'QM' initial sampling)")
parser.add_argument("--TVIBb", type=float, help="Vibrational temperature of atoms in second group (supply if using the 'thermal' initial sampling)")
parser.add_argument("--TROTb", type=float, help="Rotational temperature of atoms in second group (supply if using the 'thermal' initial sampling)")
args = vars(parser.parse_args())

########################################################################################

# Get the various arguments
Qfile = args["initialGeometryFile"]
input_path = args["psi4inputFile"]
output_path = args["outputDir"]

Pfile = args["initialMomentaFile"]
atomsInFirstGroup = args["atomsInFirstGroup"]
ce = args["collisionEnergy"]
b = args["impactParameter"]
dCM = args["centerOfMassDistance"]

isotopeMassesFile = args["isotopeMassesFile"]
Nsteps = args["production"]
Nprint = args["interval"]
dt = args["time_step"]

n_threads = args["n_threads"]
if (n_threads is None): n_threads = 1

samplingMethod      = [ args["INITQPa"], args["INITQPb"] ]
vibrationSampling = []
rotationSampling = []
if (samplingMethod[0] == "semiclassical"):
  vibrationSampling.append(args["NVIBa"])
  rotationSampling.append(args["NROTa"])
else:
  vibrationSampling.append(args["TVIBa"])
  rotationSampling.append(args["TROTa"])
if (samplingMethod[1] == "semiclassical"):
  vibrationSampling.append(args["NVIBb"])
  rotationSampling.append(args["NROTb"])
else:
  vibrationSampling.append(args["TVIBb"])
  rotationSampling.append(args["TROTb"])

# Note:
# Right now, all arguments are mandatory (even though this does
# not raise a warning) except for isotopeMassesFile

# Adjust the maximum interatomic distance allowed
# for the simulation
r2threshold = 900.0
if ((b is not None) and (dCM is not None) and (1.2*(b**2 + dCM**2) > r2threshold)):
  r2threshold = 1.2*(b**2 + dCM**2)

########################################################################################

# Get the model ready
calc = psi4calculator(input_path,n_threads=n_threads)

# To conform to VENUS, we are going to keep the units
# in kcal/mol and Angstroms (which the model was
# originally trained on)
calc.E_to_eV = units.Ha
calc.Ang_to_R = units.Ang
calc.F_to_eV_Ang = (units.Ha / units.Bohr)

# Read in the geometry; set it in the "calculator"
mol = read(Qfile)
mol.set_calculator(calc)
mol.calc = calc
mol._calc = calc

# Get the output trajectory file ready
trajfile = os.path.join(output_path, "production.traj")

# If the masses are given, update the masses
# Note: this must be done BEFORE setting the momenta
if not (isotopeMassesFile is None):
  massFile = open(isotopeMassesFile,"r")
  newMasses = massFile.readlines()
  massFile.close()
  if (len(newMasses) != len(mol)): #.masses)):
    raise ValueError("Number of masses provided in the isotope mass file does not match the input XYZ")
  mol.set_masses([float(i) for i in newMasses])

# If a momenta file is given, read in the momenta
# Read it in as a geometry, and then set it into the molecule
if not (Pfile is None):
  frame = read(Pfile)
  p = [atom.position for atom in frame]
  mol.set_momenta(p)

  masses=mol.get_masses()
  v = [p[i]/masses[i] for i in range(len(p))]
  mol.set_velocities(v)

# If there is no momenta file, then try and make momenta
# by assuming a bimolecular reaction
else:
    if ((atomsInFirstGroup is None) or (ce is None) or (b is None) or (dCM is None)):
      raise ValueError("No initial momenta conditions or file are provided")

    atomsInFirstGroup = [int(i)-1 for i in atomsInFirstGroup.split()]

    Natoms = len(mol)
    atomsInSecondGroup = []
    for i in range(Natoms):
      if (i not in atomsInFirstGroup): atomsInSecondGroup.append(i)

    print("")
    print("GEOMETRY INPUT")
    print("  Input geometry file: ", Qfile)
    print("Atoms in  first group: ", atomsInFirstGroup)
    print("Atoms in second group: ", atomsInSecondGroup)
    print("SAMPLING INPUT")
    print("  Input momenta file: ", Pfile)
    if (Pfile is None):
      print("     Group A sampling method: ", samplingMethod[0])
      print("     Group A       vibration: ", vibrationSampling[0])
      print("     Group A        rotation: ", rotationSampling[0])
      print("     Group B sampling method: ", samplingMethod[1])
      print("     Group B       vibration: ", vibrationSampling[1])
      print("     Group B        rotation: ", rotationSampling[1])
      print("        Impact parameter (A): ", b)
      print("      Initial separation (A): ", dCM)
      print("  Collsion energy (kcal/mol): ", ce)
    else:
      print("  Input momenta file: ", Pfile)

    print("##############################################################")
    print("")

    # Sample the internal positions and momenta of each of
    # the two molecules
    sampler = initialSampling(mol,atomsInFirstGroup,optimize=True,
                      optimization_file=os.path.join(output_path, "optimization.traj"),
                      samplingMethodA=samplingMethod[0],vibrationalSampleA=vibrationSampling[0],rotationalSampleA=rotationSampling[0],
                      samplingMethodB=samplingMethod[1],vibrationalSampleB=vibrationSampling[1],rotationalSampleB=rotationSampling[1])

    print("Sampling internal degrees of freedom...")
    sampler.sampleRelativeQP()

    print("Sampling relative degrees of freedom...")
    sampler.sampleAbsoluteQP(ce,dCM=dCM,b=b)

########################################################################################


# Run MD with constant energy using the velocity verlet algorithm
dyn = VelocityVerlet(mol, dt * units.fs, trajectory=trajfile) 

# A function to print the potential, kinetic and total energy
def printenergy(a):
    epot = a.get_potential_energy() / (units.kcal/units.mol)
    ekin = a.get_kinetic_energy() / (units.kcal/units.mol)
    print('@Epot = %.3f  Ekin = %.3f (T=%3.0fK)  '
          'Etot = %.3f  kcal/mol' % (epot, ekin, ekin / (len(a) * 1.5 * 8.617281e-5), epot + ekin))

# A function to see if any interatomic distance is > 20 
def checkGeneralReactionProgress(a):
    Natoms = len(a)
    stop_flag = False
    for i1 in range(1,Natoms):
      for i2 in range(0,Natoms-1):

        # Check the squared interatomic distances if any
        # are greater than r2threshold (default is 20^2 = 400)
        r = sum([(a.positions[i1][i] -
                  a.positions[i2][i])**2 for i in range(3)])
        if (r > r2threshold):
          stop_flag = True
          break

    return stop_flag

# Now run the dynamics
printenergy(mol)
for i in range(Nsteps):
  dyn.run(Nprint)
  printenergy(mol)
  stop_flag = checkGeneralReactionProgress(mol)
  if (stop_flag): break

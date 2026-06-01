import numpy as np
from ase.io import read, write
from ase import units
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from ase.calculators.calculator import Calculator, all_changes
from ase import build
import time
import torch
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


#######################################################################
# --- USER SETTINGS ---
#######################################################################

# Input structures (same number of atoms, same ordering)
xyz_A = "stateA.xyz"
xyz_B = "stateB.xyz"

# Foundational MLFF calculators for state A and B
# Example:
# from ase.calculators.mace import MACECalculator
# MLFF_A = lambda: MACECalculator(model_paths=["model_A.pt"])
# MLFF_B = lambda: MACECalculator(model_paths=["model_B.pt"])
#
# Replace these with your foundational MLFFs:

from mace.calculators import MACECalculator
MLFF_A = lambda: MACECalculator(model_paths='./mace-mpa-0-medium.model',  enable_cueq=True, device='cuda')
MLFF_B = lambda: MACECalculator(model_paths='./mace-mpa-0-medium.model',  enable_cueq=True, device='cuda')

temperature = 300                 # Kelvin
dt = 0.5 * units.fs                 # MD timestep
total_steps = 20000               # total simulation
equilibration_steps = 10000        # discard initial samples
sample_interval = 10              # how often to record dU/dλ

lambda_values = np.linspace(0, 1, 3)   # 0.0, 0.1, ..., 1.0


#######################################################################
# --- ALCHEMICAL CALCULATOR ---
#######################################################################

class AlchemicalMLFF(Calculator):
    """
    Generic alchemical calculator that mixes two arbitrary ASE MLFFs.
    Uses linear interpolation of energies & forces:
        U(λ) = (1−λ) U_A + λ U_B
    """
    implemented_properties = ["energy", "forces"]
    nolabel = True

    def __init__(self, MLFF_A_factory, MLFF_B_factory, lam, **kwargs):
        super().__init__(**kwargs)
        self.lam = lam
        self.calcA = MLFF_A_factory()
        self.calcB = MLFF_B_factory()

    def calculate(self, atoms=None, properties=["energy", "forces"], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        # Energies
        eA = self.calcA.get_potential_energy(atoms)
        eB = self.calcB.get_potential_energy(atoms)

        # Forces
        fA = self.calcA.get_forces(atoms)
        fB = self.calcB.get_forces(atoms)

        lam = self.lam

        # Interpolated results
        self.results["energy"] = (1 - lam) * eA + lam * eB
        self.results["forces"] = (1 - lam) * fA + lam * fB

        # Tight TI quantity: ∂U/∂λ = U_B − U_A
        self.results["dU_dlambda"] = eB - eA


#######################################################################
# --- LOAD STRUCTURES AND INTERPOLATE GEOMETRY ---
#######################################################################

atoms_A = read(xyz_A)
atoms_A.set_cell([45.329860, 21.290961, 9.46])
atoms_A.set_pbc(pbc=True)
MaxwellBoltzmannDistribution(atoms_A, temperature)
atoms_B = read(xyz_B)
atoms_B.set_cell([45.329860, 21.290961, 9.46])
atoms_B.set_pbc(pbc=True)
MaxwellBoltzmannDistribution(atoms_B, temperature)

if len(atoms_A) != len(atoms_B):
    raise ValueError("State A and B must have the same number of atoms. "
                     "If you need atom mapping or padding, tell me.")

# Simple geometric interpolation (species come from stateA)
def interpolate_geometry(lam):
    atoms = atoms_A.copy()
    atoms.positions = (1 - lam) * atoms_A.positions + lam * atoms_B.positions
    return atoms


#######################################################################
# --- RUN MD AT EACH λ ---
#######################################################################

dU_dlambda_avg = {}

for lam in lambda_values:
    print(f"\n=== Running λ = {lam:.2f} ===")

    atoms = interpolate_geometry(lam)

    atoms.calc = AlchemicalMLFF(MLFF_A, MLFF_B, lam)

    dyn = NoseHooverChainNVT(atoms, dt, temperature, taut=100 * dt)

    samples = []

    def record(a=atoms):
        samples.append(a.calc.results["dU_dlambda"])

    dyn.attach(record, interval=sample_interval)
    dyn.run(total_steps)

    # Average over production region
    dU_dlambda_avg[lam] = np.mean(samples[equilibration_steps:])
    print(f"⟨dU/dλ⟩ = {dU_dlambda_avg[lam]:.6f} eV")


#######################################################################
# --- THERMODYNAMIC INTEGRATION ---
#######################################################################

lams = np.array(list(dU_dlambda_avg.keys()))
vals = np.array(list(dU_dlambda_avg.values()))
lam_points = sorted(dU_dlambda_avg.keys())
if len(lam_points)!=3:
    raise ValueError ("Simpson's 1/3 rule requires exactly 3 lambda points.")
lam0, lam1, lam2 =lam_points
f0=dU_dlambda_avg[lam0]
f1=dU_dlambda_avg[lam1]
f2=dU_dlambda_avg[lam2]
h=lam2-lam0
dFT=np.trapezoid(vals, lams)*96.49
dFS = (h/6.0)*(f0+4*f1+f2)*96.49

print("\n===============================================")
print(f"ΔF(A → B) = {dFS:.2f} kJ/mol, trapezoid= {dFT:.2f} kJ/mol.")
print("===============================================")
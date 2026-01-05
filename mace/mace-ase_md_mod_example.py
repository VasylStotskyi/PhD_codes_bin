from ase import units
from ase.md.nose_hoover_chain import NoseHooverChainNVT
from ase.io import read, write
from ase.io import Trajectory
import numpy as np
import time
import torch
from mace.calculators import MACECalculator
from ase import build
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

energy_list=[]
calculator = MACECalculator(model_path='./mace-mpa-0-medium.model')
# return the default medium ASE calculator equivalent to mace_mp(model="medium")
#macemp = mace_mp(model="large") # return a larger model
#macemp = mace_mp(model="https://tinyurl.com/y7uhwpje") # downlaod the model at the given url
#mace
init_conf = read('frame16000_17.xyz', '0')
init_conf.set_cell([45.3437, 21.2845, 9.46])
init_conf.set_pbc(pbc=True)
init_conf.set_calculator(calculator)
# Initialize velocities.
temperature_K = 300 # Initial temperature in K
MaxwellBoltzmannDistribution(init_conf, temperature_K * units.kB)

# Set up the Nose-Hoover dynamics engine for NVT ensemble.
dyn = NoseHooverChainNVT(init_conf, 0.5 * units.fs, temperature_K=300, tdamp=100*units.fs)
def write_frame():
        kinetic_energy=dyn.atoms.get_kinetic_energy()
        temperature=dyn.atoms.get_temperature()
        potential_energy=dyn.atoms.get_potential_energy()
        total_energy=dyn.atoms.get_total_energy()
        dyn.atoms.write('md_300K_mace.xyz', format='xyz', append=True)
        energy_list.append([kinetic_energy, temperature, potential_energy,total_energy])
#traj=Trajectory('md_300K_mace.xyz','w', format='xyz', append=True)
dyn.attach(write_frame, interval=2)
n_steps = 10 # Number of steps to run
dyn.run(n_steps)
ener=open('md_300K_mace.ener', "w")
for line in energy_list:
        print('\t'.join(str(x) for x in line), file=ener)
ener.close()


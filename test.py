import os
import numpy as np
from ase.io import read
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.calculators.lj import LennardJones
ASE_VASP_COMMAND="mpirun -np "+str(16)+" vasp_std"

os.environ["VASP_PP_PATH"]=os.path.join(os.getcwd(),"VASP_PP")

calc = Vasp(xc='PBE',
            encut=400,
            ediff=1e-4,
            ismear=0,
            sigma=0.02,
            kspacing=0.2,
            gamma=True,
            nelm=200,
            restart=None,            
            command=ASE_VASP_COMMAND,
            directory="a/vasp_run",
            txt=None)

atoms=read("B3_sub_minimal/traj_seed_0_0.res")
atoms.set_calculator(calc)
print(atoms.get_potential_energy())

positions=atoms.get_positions()
cell=[[2.087777, 0.0, 0.0], [-5.632629886459229e-07, 2.933871999999946, 0.0], [-1.0438533985624054, -1.4669148420831921, 8.632346199162008]]
atoms=Atoms("B"*len(atoms),positions=positions,cell=cell)
atoms.pbc=True
atoms.set_calculator(calc)
print(atoms.get_potential_energy())
print(type(atoms.get_positions()))
print(type(np.array(atoms.get_cell())))
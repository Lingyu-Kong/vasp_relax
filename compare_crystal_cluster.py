import os
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

## Treat as a Crystal
atoms.set_calculator(calc)
energy=atoms.get_potential_energy()
print("Treatment as a Crystal: ",energy)

## Treat as a Cluster
positions=atoms.get_positions()
pos_mean=positions.mean(axis=0)
positions=positions-pos_mean
cluster=Atoms("B"+str(len(positions)),positions=positions)
calc=LennardJones(rc=500)
cluster.set_calculator(calc)
energy=cluster.get_potential_energy()
print("Treatment as a Cluster: ",energy)

os.system("rm -r a/vasp_run")
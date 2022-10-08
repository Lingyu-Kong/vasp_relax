from ase.io.trajectory import Trajectory
from ase.io import read,write

traj = Trajectory('B3_sub_minimal[0, 2].traj')
print(len(traj))
atoms = traj[-1]
print(atoms.get_potential_energy())
print(atoms.get_forces())
print(atoms.get_stress())
print(atoms.get_positions())
print(atoms.get_cell())

atoms = traj[-2]
print(atoms.get_potential_energy())
print(atoms.get_forces())
print(atoms.get_stress())
print(atoms.get_positions())
print(atoms.get_cell())
    
# atoms=read('./B8_minimal/traj_seed_0.res')
# print(atoms.get_scaled_positions())

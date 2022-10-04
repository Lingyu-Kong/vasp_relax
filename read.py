from ase.io.trajectory import Trajectory

traj = Trajectory('./B8_minimal_shake/traj_seed_0.traj')
print(len(traj))
atoms = traj[-1]
print(atoms.get_potential_energy())
print(atoms.get_forces())
print(atoms.get_stress())
print(atoms.get_positions())
print(atoms.get_cell())

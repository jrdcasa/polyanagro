import polyanagro as pag
import numpy as np
from ext_libc.c_unit_bond_vectors import unit_bond_vectors

# Logger
log = pag.init_logger("Output", fileoutput="example_unitvectors.log",
                          append=False, inscreen=False)
# End Logger
# Setup Trajectory to analyze
xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
filename_tpr = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
filename_psf = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/00-GENERATE/namd_out.psf"
stride = 1
trj_small = pag.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=log)
trj_small.topology.assign_backbone_atoms_psf()
# End Setup Trajectory to analyze

natoms = trj_small.topology._natoms
nmols_array, neigbors = trj_small.topology.get_array_mols_neigh()
nchains = len(nmols_array)
nframes = trj_small.get_numframes()
iatch = trj_small.topology._iatch

# Get backbone atoms
all_bb_bonds, nbonds_bb_perch = trj_small.topology.get_all_bb_bonds()
all_bb_bonds = np.array(all_bb_bonds, dtype=np.int32)
coords_t0_wrapped = trj_small.universe.trajectory[0].positions
box_dimensions = trj_small.universe.trajectory[0].dimensions[0:3]
nbonds_bb_max = max(nbonds_bb_perch.values())
coords_unwrap = pag.unwrap(coords_t0_wrapped, nmols_array, neigbors,
                           box_dimensions, iframe=0)

# Preparing unit bond vector
uux = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
uuy = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
uuz = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)

unit_bond_vectors(nchains, all_bb_bonds, coords_unwrap, iatch, uux, uuy, uuz)

print(uux)
print("======")
print(uuy)
print("======")
print(uuz)
# print(uuz)
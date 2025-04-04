#from polyanagro.Topology import Topology
from polyanagro.ExtTrajectory import Trajectory
#from polyanagro.utils.Logger import LogCreate
#from polyanagro.Chain_Statistics import Chain_Statistics
#from polyanagro.Calculations import unwrap_purepython, unwrap
#import numpy as np
#from polyanagro.utils import total_size
#import datetime

# Create the logger
LogCreate(["console"])

# Input files
xtc1 = "../../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
xtc2 = "../../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
xtc3 = "../../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
filename_tpr = "../../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
stride = 1

# Setup Trajectory
trj = Trajectory([xtc1, xtc2, xtc3])

# Setup topology
top = Topology()
top.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
trj.set_top(top)

# Create a list with atoms to calculate the end-to-end distances
end2end_atoms = []
nmols_array, neigbors = top.get_array_mols_neigh()
nchains = len(nmols_array)
listend2end = []
natch = 24
natbb = 20
for ich in range(nchains):
    ihead = 0 + ich*natch
    itail = (natbb-1) + (ich*natch)
    listend2end.append([ich, ihead, itail])

# Calculations
objcalc = Chain_Statistics(trj, top, dt=20, stride=stride)
#objcalc.calculate(diroutput="./", listendtoend=listend2end, acfE2E=True, distributions=True)
objcalc.calculate(diroutput="./", listendtoend=listend2end, acfE2E=False, distributions=True)
#objcalc.end2end(filename="Ree.dat")

# print(total_size(c))
# print("%d bytes"%(c.size*c.itemsize))
# pass
from polyanagro.Topology import Topology
from polyanagro.ExtTrajectory import Trajectory
from polyanagro.Logger import LogCreate
from polyanagro.Calculations import Calculations

# Create the logger
LogCreate(["console"])

# Input files
xtc1 = "/media/jramos/JRD_05/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
xtc2 = "/media/jramos/JRD_05/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
xtc3 = "/media/jramos/JRD_05/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
filename_tpr = "/media/jramos/JRD_05/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
stride = 1000

# Setup Trajectory
trj = Trajectory([xtc1, xtc2, xtc3])

# Setup topology
top = Topology()
top.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
trj.set_top(top)

# Calculations
objcalc = Calculations(trj, top, dt=20, stride=stride)
objcalc.radius_of_gyration()
#objcalc.end2end(filename="Ree.dat")

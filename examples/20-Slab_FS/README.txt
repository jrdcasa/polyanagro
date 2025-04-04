gmx trjconv -f ORIGINALS/traj.part0001_0-100ns.trr  -s new_topo.tpr -o traj_cut.trr -skip 10

Command line:
  gmx check -f traj_cut.trr

Checking file traj_cut.trr
trr version: GMX_trn_file (single precision)
Reading frame       0 time    0.000   
# Atoms  73440
Last frame         50 time 100000.000   


Item        #frames Timestep (ps)
Step            51    2000
Time            51    2000
Lambda          51    2000
Coords          51    2000
Velocities      51    2000
Forces           0
Box             51    2000



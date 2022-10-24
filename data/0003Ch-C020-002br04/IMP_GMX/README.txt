module load GROMACS
gmx trjcat -f ../RUN-001/traj_comp.xtc ../RUN-002/traj_comp.part0002.xtc ../RUN-003/traj_comp.part0003.xtc -o trj_wrap_full.xtc
gmx trjconv -f trj_wrap_full.xtc -s ../RUN-001/topol.tpr -pbc whole -o trj_unwrap_full.xtc
vmd ../RUN-001/confout.gro trj_unwrap_full.xtc

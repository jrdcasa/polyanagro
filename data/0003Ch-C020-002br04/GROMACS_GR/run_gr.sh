
source /opt/gromacs_2021.3/nogpu/bin/GMXRC.bash

gmx_cpu rdf -f ../RUN-001/traj_comp.xtc  -s ../RUN-001/topol.tpr -n index_all.ndx -ref setA -sel setB -o rdf_all.xvg

gmx_cpu rdf -f ../RUN-001/traj_comp.xtc  -s ../RUN-001/topol.tpr -n index_bychain.ndx -ref chain1 -sel chain2 -o rdf_chain1_chain2.xvg

rm \#*



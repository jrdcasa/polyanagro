exec(open('/usr/local/Modules/init/python.py').read())
module('list')
module('available')
module('load', 'gromacs-2018.3')

import gromacs
gromacs.release()
print(gromacs.release())



# Input files
xtc1 = "../../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
xtc2 = "../../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
xtc3 = "../../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
filename_tpr = "../../data/0003Ch-C020-002br04/RUN-001/topol.tpr"

gromacs.g_dist(s=filename_tpr, f=xtc1, oh='dist.xvg',oall="alldist.xvg", input="0")
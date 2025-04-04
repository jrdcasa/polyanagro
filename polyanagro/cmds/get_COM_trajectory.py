import MDAnalysis as mda
import numpy as np
import datetime
import argparse
import os
import utils
import sys
import topology


# =============================================================================
def parse_arguments():

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    desc = """Transform an atomic trajectory to a COM trajectory.
    This is part of the polyanagro library"""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.\n"
                             "Allowed trajectories are XTC and DCD.",
                        action="store", required=True, default=None)

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="info_comtrj.log")

    parser.add_argument("--topo", dest="topo",
                        help="A topology file in tpr, data or pdb format.\n"
                             "tpr --> GROMACS, dat --> LAMMPS, pdb --> OTHERS",
                        action="store", metavar="TPR|DATA|PDB", required=True)

    parser.add_argument("--input_nojump", dest="input_nojump",
                        help="If True, the input trajectory is already unwrapped without jumps, otherwise"
                             "the program unwrap without jumps the input trajectory.",
                        action="store_true", required=False)

    args = parser.parse_args()

    for itrj in args.traj:
        if not os.path.isfile(itrj):
            print("\nERROR: File {} does not exist\n".format(itrj))
            exit()

    if not os.path.isfile(args.topo):
        print("\nERROR: File {} does not exist\n".format(args.topo))
        exit()

    return args

# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                         Center mol in a  Trajectory 
              ----------------------------------------------

                                Version {}

                              Dr. Javier Ramos
                      Macromolecular Physics Department
                Instituto de Estructura de la Materia (IEM-CSIC)
                               Madrid (Spain)

        This utility is part of the polyanagro library. Polyanagro is an 
        open-source python library to analyze simulations of polymer systems.

        This software is distributed under the terms of the
        GNU General Public License v3.0 (GNU GPLv3). A copy of
        the license (LICENSE.txt) is included with this distribution.

    ***********************************************************************
        """.format(version)

    print(msg) if logger_log is None else logger_log.info(msg)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

    m = "\t\tStart Job at {} ============".format(now)
    print(m) if logger_log is None else logger_log.info(m)

    m1 = ""
    for item in sys.argv[1:]:
        m1 += " {}".format(item)
    m = "\n\t\tCommand line: \n"
    m += "\t\t\tpython {}".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    m += "\t\t\t         or\n"
    m += "\t\t\tcenter_trajectory".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    print(m) if logger_log is None else logger_log.info(m)


# =============================================================================
def main_app():

    import polyanagro as pag

    # Parse arguments
    args = parse_arguments()
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)
    # Write header and arguments
    print_header(pag.version.__version__, log)

    # Load trajectory
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    trj = topology.ExtTrajectory(args.traj, topfile=args.topo, logger=log)

    # Print Info
    natoms = trj.universe.trajectory.n_atoms
    nmols = len(trj.topology._nmols)
    nframes = trj.universe.trajectory.n_frames
    m = "\t**** TRAJECTORY INFO ****\n"
    m += "\t  Number of trj files         : {}\n".format(len(args.traj))
    m += "\t  Number of atoms             : {}\n".format(natoms)
    m += "\t  Number of molecules (chains): {}\n".format(nmols)
    m += "\t  Number of atoms/molecule    : {}\n".format(natoms/nmols)
    m += "\t  Number of frames            : {}\n".format(nframes)
    m += "\t  Initial time (ps)           : {}\n".format(trj.universe.trajectory.time)
    m += "\t  Time step (ps)              : {}\n".format(round(trj.universe.trajectory.dt))
    m += "\t  Total time (ps)             : {}\n".format(round(trj.universe.trajectory.dt*
                                                              (nframes-1)))
    m += "\t**** End TRAJECTORY INFO ****"
    print(m) if log is None else log.info(m)

    #For each frame in the trajectory
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\t Calculating center of mass ... ({})\n".format(now)
    print(m) if log is None else log.info(m)

     # No jump trajectory
    dimensions_trj = []
    dt_trj = []
    for ts in trj.universe.trajectory:
        dt_trj.append(ts.dt)
        dimensions_trj.append(ts.dimensions)

    if args.input_nojump:
        trj.write_trajectory("tmp_nojump", nojump=False, format_trj="xtc")
    else:
        trj.write_trajectory("tmp_nojump", nojump=True, format_trj="xtc")
    trj = topology.ExtTrajectory("tmp_nojump.xtc", topfile=args.topo, logger=log)
    with open("com_trajectory.xyz", "w") as fxyz:
        iframe = 0
        indx_iframe = 0
        start_time = datetime.datetime.now()
        for ts in trj.universe.trajectory:
            if indx_iframe % 100 == 0:
                msg = "\tUnwrap Nojump Frame {0:9d} of {1:9d}".format(indx_iframe, nframes)
                mid_time = datetime.datetime.now()
                elapsed_time = mid_time - start_time
                msg += "\ttime: {0:s} seconds".format(str(elapsed_time.total_seconds()))
                print(msg) if log is None else log.info(msg)
            fxyz.writelines("{}\n".format(nmols))
            fxyz.writelines("Frame {}\n".format(iframe))
            #For each molecule (chain) in the frame
            for imol_list in trj.topology._nmols:
                # Select atoms for imol
                select_str = "index " + " ".join(map(str, imol_list))
                atg = trj.universe.select_atoms(select_str)
                com = atg.center_of_mass()
                fxyz.writelines("C {} {} {}\n".format(com[0], com[1], com[2]))
            indx_iframe += 1

    # Load XYZ trajectory
    u = mda.Universe("com_trajectory.xyz", "com_trajectory.xyz", format="XYZ")

    # Write to XTC
    iframe = 0
    box_dims = [ 52.4992,   52.4992,   52.4992, 90.0, 90.0, 90.0]
    with mda.Writer("com_trajectory.xtc", n_atoms=u.atoms.n_atoms) as w:
        for ts in u.trajectory:
            ts.dimensions = dimensions_trj[iframe]
            ts.time = dt_trj[iframe]
            w.write(u.atoms)
            iframe += 1
    # Set the first frame
    first_frame = u.trajectory[0]  # Select the first frame
    first_frame.dimensions = dimensions_trj[0]  # Assign box dimensions
    first_frame.time = dt_trj[0]

    # Write the first frame to GRO format
    with mda.Writer("com_trajectory.gro", n_atoms=u.atoms.n_atoms) as w:
        w.write(u.atoms)  # Write AtomGroup

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)

# =============================================================================
if __name__ == "__main__":

    main_app()


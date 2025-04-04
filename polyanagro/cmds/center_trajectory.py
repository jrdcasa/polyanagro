import MDAnalysis as mda
import MDAnalysis.transformations as trans
import datetime
import argparse
import os
import utils
import sys


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

    desc = """Get information from a MD trajectory.
    This is part of the polyanagro library"""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.\n"
                             "Allowed trajectories are XTC and DCD.",
                        action="store", required=True, default=None)

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="info_centertrj.log")

    parser.add_argument("--topo", dest="topo",
                        help="A topology file in tpr, data or pdb format.\n"
                             "tpr --> GROMACS, dat --> LAMMPS, pdb --> OTHERS",
                        action="store", metavar="TPR|DATA|PDB")

    parser.add_argument("--sel", dest="selection",
                        help="Selection to center in the simulation box",
                        action="store", required=True)

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
    m = "\t\t Loading trajectory ({})\n".format(now)
    for item in args.traj:
        m += "\t\t\t {}\n".format(item)
    print(m) if log is None else log.info(m)
    u = mda.Universe(args.topo, args.traj, in_memory=True)

    center_mol = u.select_atoms(args.selection)
    all_atoms = u.select_atoms("all")
    not_selection = u.select_atoms('not '+args.selection)
    not_water = u.select_atoms("not resname SOL")

    # Print Info
    m = "\t**** TRAJECTORY INFO ****\n"
    m += "\t  Number of trj files        : {}\n".format(len(args.traj))
    m += "\t  Number of atoms            : {}\n".format(u.trajectory.n_atoms)
    m += "\t  Number of frames           : {}\n".format(u.trajectory.n_frames)
    m += "\t  Time step (ps)             : {}\n".format(round(u.trajectory.dt))
    m += "\t  Total time (ps)            : {}\n".format(round(u.trajectory.dt*(u.trajectory.n_frames-1)))
    m += "\t**** End TRAJECTORY INFO ****"
    print(m) if log is None else log.info(m)

    # Center protein in box
    # https://userguide.mdanalysis.org/stable/examples/transformations/center_protein_in_box.html
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t**** Unwrapping the selection **** ({})\n".format(now)
    print(m) if log is None else log.info(m)
    transforms = [trans.unwrap(center_mol),
                  trans.center_in_box(center_mol, wrap=True, center='geometry'),
                  trans.wrap(not_selection)]

    u.trajectory.add_transformations(*transforms)

    # Write new trajectory
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t**** Writing trajectory **** ({})\n".format(now)
    print(m) if log is None else log.info(m)
    with mda.Writer("traj2-Achain.xtc", all_atoms.n_atoms) as W:
        for ts in u.trajectory:
            W.write(all_atoms)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)

# =============================================================================
if __name__ == "__main__":

    main_app()


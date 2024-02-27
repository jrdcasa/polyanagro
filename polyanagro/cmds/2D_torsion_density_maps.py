import utils
import os
import sys
import argparse
import datetime
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

    desc = """Bonded distributions for atomistic and coarse-grained models."""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.",
                        action="store", required=True, default=None)

    parser.add_argument("--topo", dest="topo",
                        help="A topology file in tpr, data or pdb format.\n"
                             "tpr --> GROMACS, dat --> LAMMPS, pdb --> OTHERS",
                        action="store", required=True, metavar="TPR|DATA|PDB")

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="pol_2DMap.log")

    parser.add_argument("--phipsi", dest="phipsi", nargs=2,
                        help="""A list with two labels contained in the file dihedrals_data_dist.ndx.""",
                        action="store", required=True, default=None)

    parser.add_argument("--stride", dest="stride",
                        help="Take a frame each stride frames, for example 10",
                        action="store", required=False, default=1)

    parser.add_argument("--unwrap", dest="isunwrap", type=str2bool,
                        help="If True the coordinates provided are unwrapped",
                        required=True, metavar="True or False, 1 or 0")

    parser.add_argument("--half_dihedral", dest="half_dihedral",
                        help="If True the range of dihedrals is [-180..180] "
                             "otherwise is [0..360]",
                        required=False, action="store_true")

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
                  Generate a 2DMap for a pair of dyhedral types
        -------------------------------------------------------------

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
    m += "\t\t\t{}".format(os.path.split(sys.argv[0])[1])
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
    trj = topology.ExtTrajectory(args.traj, topfile=args.topo, logger=log)

    # Check if dihedral_data_dist.ndx exists
    ndxfilename = "dihedral_data_dist.ndx"
    if not os.path.isfile(ndxfilename):
        m = "\n\t\t The file dihedral_data_dist.ndx cannot be found.\n"
        m += "\t\t You can produce this file with the bond_distribution generate tool.\n"
        m += "\t\t presents in this polyanagro suite.\n"
        print(m) if log is None else log.error(m)
        exit()

    pairdihlist = args.phipsi
    isunwrap = args.isunwrap

    objmap = pag.Map_Density(trj, pairdihlist, log=log)
    objmap.calculate(begin=0, unwrap_pbc=isunwrap,
                     ndx_filename = ndxfilename, half_dihedral=args.half_dihedral)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)


# =============================================================================
if __name__ == "__main__":

    main_app()
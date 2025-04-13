import utils
import os
import sys
import argparse
import datetime
import topology
import time

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

    desc = """Mean Squared displacement for atomistic and coarse-grained models."""

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
                        action="store", required=False, default="MSD.log")

    parser.add_argument("--nojump", dest="nojump", type=str2bool,
                        help="If True the coordinates provided are unwrapped without jumps",
                        required=True, metavar="True or False, 1 or 0")

    parser.add_argument("--com", dest="com", type=str2bool,
                        help="If True the g3(t) function is calculated",
                        required=True, metavar="True or False, 1 or 0")

    parser.add_argument("--start", dest="start", type=int,
                        help="The starting frame for the trajectory.",
                        action="store", required=False, default=0)

    parser.add_argument("--end", dest="end", type=int,
                        help="The stopping frame for the trajectory.",
                        action="store", required=False, default=-1)

    parser.add_argument("--stride", dest="stride", type=int,
                        help="The step size for the trajectory.",
                        action="store", required=False, default=1)

    parser.add_argument("-o", "--output", dest="output",
                        help="Name of the output for the MSD data.",
                        action="store", required=False, default="msd.dat")

    parser.add_argument("--method", dest="method",
                        help="The method to calculate the MSD. "
                             "Values: classic_opt_cython, classic_multitau,  msd_fftw3_fast, msd_fftw3_cython."
                             "Default: msd_fftw3_cython",
                        action="store", required=False, default="msd_fftw3_cython")

    parser.add_argument("--startingfactor", dest="startingfactor", type=int,
                        help="Number of starting points of 'independent' time series",
                        action="store", required=False, default=None)

    args = parser.parse_args()

    for itrj in args.traj:
        if not os.path.isfile(itrj):
            print("\nERROR: File {} does not exist\n".format(itrj))
            exit()

    if not os.path.isfile(args.topo):
        print("\nERROR: File {} does not exist\n".format(args.topo))
        exit()

    methods_available = ["classic_opt_cython", "classic_multitau",  "msd_fftw3_fast", "msd_fftw3_cython"]
    if args.method not in methods_available:
        print("\nERROR: Method is not implemented\n")
        print("Methods available are: ")
        for item in methods_available:
            print(" {}".format(item))
        print("\n")
        exit()

    if args.method == "classic_multitau":
        if args.startingfactor is None:
            print("\nERROR: classic_multitau requires a value for --startingfactor\n")
            print("\n")
            exit()

    return args


# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                     Calculate mean square displacement 
                       using different implementations
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

    start_time = time.time()

    # Parse arguments
    args = parse_arguments()
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)
    # Write header and arguments
    print_header(pag.version.__version__, log)
    # Load trajectory
    trj = topology.ExtTrajectory(args.traj, topfile=args.topo, logger=log)

    if not args.nojump:
        msg = "\t\t The original trajectory is not unwrapped and/or jumps are not corrected.\n"
        msg += "\t\t Nojump option: {}. Provided trajectory is wrapped.\n".format(args.nojump)
        msg += "\t\t  \n"
        print(msg) if log is None else log.info(msg)
        filename_nojump = os.path.splitext(os.path.split(args.traj[0])[-1])[0]+"_nojump.xtc"
        trj.write_trajectory(filename_nojump, pbc=True, nojump=True, format_trj=None)
        trj = topology.ExtTrajectory(filename_nojump, topfile=args.topo, logger=log)
        args.nojump = 1

    # objmsd = pag.MSD(trj, args.nojump, method="classic_py", outputname="msd_classic.dat", start=args.start,
    #                 end=args.end, step=args.stride, logger=log)
    #
    if args.method == "classic_opt_cython":
        objmsd = pag.MSD(trj, args.nojump, method="classic_opt_cython", outputname=args.output, start=args.start,
                         end=args.end, step=args.stride, com=args.com, logger=log)
    elif args.method == "classic_multitau":
        objmsd = pag.MSD(trj, args.nojump, method="classic_multitau", outputname=args.output, start=args.start,
                         end=args.end, step=args.stride, com=args.com, startingfactor=args.startingfactor, logger=log)
    elif args.method == "msd_fftw3_fast":
        objmsd = pag.MSD(trj, args.nojump, method="msd_fftw3_fast", outputname=args.output, start=args.start,
                         end=args.end, com=args.com, step=args.stride, logger=log)
    elif args.method == "msd_fftw3_cython":
        objmsd = pag.MSD(trj, args.nojump, method="msd_fftw3_cython", outputname=args.output, start=args.start,
                         end=args.end, com=args.com, step=args.stride, logger=log)
    else:
        msg = "\t\t The original trajectory is not unwrapped and/or jumps are not corrected.\n"
        msg += "\t\t Nojump option: {}. Provided trajectory is wrapped.\n".format(args.nojump)
        msg += "\t\t  \n"
        print(msg) if log is None else log.info(msg)

    #objmsd = pag.MSD(trj, args.nojump, method="multiple_window", logger=log)

    # Calculate duration
    end_time = time.time()
    duration = end_time - start_time
    msg = "\n\t\tExecution time: {0:.2f} seconds".format(duration)
    print(msg) if log is None else log.info(msg)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    msg = "\n\t\tJob  Done at {} ============\n".format(now)
    print(msg) if log is None else log.info(msg)


# =============================================================================
if __name__ == "__main__":

    main_app()
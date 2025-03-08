import argparse
import os
import datetime
import utils
import sys
import topology
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
import warnings

# Suppress specific MDAnalysis warnings (optional)
warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
def parse_arguments():

    desc = """Calculate the pair distribution function in polymers.
    This is part of the polyanagro library.
    For more information see "Computer Simulation of Liquids, 2017, 272-274"""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.",
                        action="store", required=True, default=None)

    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument("--tpr", dest="topo",
                        help="A topology file in tpr format.",
                        action="store")
    group2.add_argument("--psf", dest="topo",
                        help="A topology file in psf format.",
                        action="store")

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="pol_rdf.log")

    parser.add_argument("--stride", dest="stride",
                        help="Take a frame each stride frames, for example 10",
                        action="store", required=False, default=1)

    group4 = parser.add_mutually_exclusive_group(required=True)
    group4.add_argument("--sets", dest="sets", nargs=2,
                        help="Set of atoms to calculate the g(r)."
                             "By default all atoms are used in both sets to calculate the g(r) function",
                        action="store")
    group4.add_argument("--index", dest="index", type=str,
                        help="A ndx file format with two set of indices",
                        action="store")

    group3 = parser.add_mutually_exclusive_group(required=True)
    group3.add_argument("--dr", dest="dr", type=float,
                        help="Bin width in angstroms for the histogram",
                        action="store", required=False)

    group3.add_argument("--nbins", dest="nbins", type=int,
                        help="Number of bins in the histogram",
                        action="store", required=False)

    parser.add_argument("--cutoff", dest="cutoff", type=float,
                        help="Calculate the RDF out to cutoff in angstroms.",
                        action="store", required=False, default=20.0)

    parser.add_argument("--excl", dest="excl", type=int,
                        help="Depth of the neighbours exluded: "
                             "--excl 1 --> Exclude the (1-2) interations"
                             "--excl 2 --> Exclude 1-3 interations "
                             "and so on",
                        action="store", required=False, default=None)

    args = parser.parse_args()

    for itrj in args.traj:
        if not os.path.isfile(itrj):
            print("\nERROR: File {} does not exist\n".format(itrj))
            exit()

    if not os.path.isfile(args.topo):
        print("\nERROR: File {} does not exist\n".format(args.topo))
        exit()

    if not args.sets is None and len(args.sets) != 2:
        print("\nERROR: Two files with index sets are needed\n".format(args.sets))
        exit()

    return args

# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                         Pair distribution calculations 
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
    m += "\t\t\tpair_distribution".format(os.path.split(sys.argv[0])[1])
    m += m1 + "\n"
    print(m) if logger_log is None else logger_log.info(m)



# =============================================================================
def main_app():

    import polyanagro as pag

    # Parse arguments
    args = parse_arguments()
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)

    # Starting time
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tStarting at {} ============\n".format(now)
    print(m) if log is None else log.info(m)

    # Write header and arguments
    print_header(pag.version.__version__, log)
    # Load trajectory
    trj = topology.ExtTrajectory(args.traj, topfile=args.topo, logger=log)

    # Create the RDF object to calculate
    if args.nbins:
        delta_r = args.cutoff / args.nbins
    else:
        delta_r = args.dr

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tStarting at {} ============\n".format(now)
    print(m) if log is None else log.info(m)


    RDFcalc = pag.RDF(trj, dt=trj.dt, stride=args.stride, excl=args.excl,
                      delta_r= delta_r, cutoff= args.cutoff, logger=log)

    if args.sets:
        RDFcalc.define_sets(setA=args.sets[0], setB=args.sets[1])
    if args.index:
        RDFcalc.define_sets_idx(idx_file=args.index)

    RDFcalc.rdf_calc_cython()

    # End time
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\n\t\t End at {} ============\n".format(now)
    m += "\t\t Job Done!!!!"
    print(m) if log is None else log.info(m)

# =============================================================================
if __name__ == "__main__":

    main_app()
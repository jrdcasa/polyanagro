import utils
import os
import sys
import argparse
import datetime
import topology
import numpy as np


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

    desc = """Calculate the cohesive energy using GROMACS."""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from GROMACS MD simulations.",
                        action="store", required=True, default=None)

    parser.add_argument("--tpr", dest="tpr",
                        help="A topology file in tpr format.\n"
                             "tpr --> GROMACS",
                        action="store", required=True, metavar="TPR")

    parser.add_argument("--topo", dest="topo",
                        help="A topology file in topo format.\n"
                             "topo --> GROMACS",
                        action="store", required=True, metavar="TOPO")


    parser.add_argument("--gmxfullpath", dest="gmxfullpath", type=str,
                        help="Path to the program to join energy files.\n"
                             "Example: For GROMACS --> /usr/bin/gmx",
                        action="store", required=True, default=None)

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="energy_cohesive.log")

    parser.add_argument("--ndx", dest="ndx",
                        help="""Filename of the ndx file for molecules.""",
                        action="store", required=False, default="index.ndx")
    parser.add_argument("--cutoffs", dest="cutoffs", nargs=3,
                        help="[<rvdw cutoff(in nm)> <rcoulomb cutoff(in nm)> <Dispersion_Correction 1 or 0>",
                        action="store", required=False, default=None)

    parser.add_argument("--fraction_trj_avg", dest="frac_avg", type=float,
                        help="""Fraction of the trajectory discarted to calculate the averages. 
                        Example: 0.25 means that the 25% first frames are discarted in the average calculation.
                                 0.00 means that the full trajectory is considered in the average""",
                        action="store", required=False, default=0.0)

    args = parser.parse_args()

    if not os.path.isfile(args.tpr):
        print("\nERROR: The tpr {} file does not exist.\n".format(args.tpr))
        exit()

    if not os.path.isfile(args.gmxfullpath):
        print("\nERROR: The path {} to gmx does not exist.\n".format(args.gmxfullpath))
        exit()

    if args.cutoffs is None:
        # rvdw, rcoulomb, dispersioncorrect
        args.cutoffs = [1.0, 1.0, 1]

    return args


# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                  Calculate the cohesive energy using GROMACS
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
    gmxexepath = args.gmxfullpath
    tpr = args.tpr
    trjlist = args.traj
    topo = args.topo
    fraction_trj_avg = args.frac_avg
    rvdw, rcoulomb, dispersioncorrect = args.cutoffs
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=True)
    # Write header and arguments
    print_header(pag.version.__version__, log)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tStarting at {} ============\n".format(now)
    print(m) if log is None else log.info(m)

    m = "\t\t Working directory: {}\n".format(os.getcwd())
    print(m) if log is None else log.info(m)

    # Load trajectory

    molecularweight = True
    trj = topology.ExtTrajectory(args.traj, topfile=tpr, logger=log)
    # Write messages
    start_time = datetime.datetime.now()
    if molecularweight:
        m = ""
        m += "\t*** Molecular weigth...\n"
        print(m) if log is None else log.info(m)
        nmols = len(trj.topology._nmols)
        mass_by_chain = np.zeros(nmols)
        for ich in range(nmols):
            m = 0.0
            for el in trj.topology._nmols[ich]:
                m += trj.topology.mass[el]
                mass_by_chain[ich] = m
        molecularweigth_avg = np.mean(mass_by_chain)
        total_mass = np.sum(mass_by_chain)
        m = "\tAverage Molecular Weigth (g/mol) : {0:.2f}\n".format(molecularweigth_avg)
        m += "\tTotal Molecular Mass (g/mol) : {0:.2f}\n".format(total_mass)
        print(m) if log is None else log.info(m)

    obj = pag.EnergyCohesive(gmxexepath, tpr, trjlist, topo,
                             fraction_trj_avg, logger=log)
    obj.create_molindex_ndx(args.ndx)

    obj.get_isolatedmol_energy(rvdw=rvdw, rcoulomb=rcoulomb, dispersioncorrect=dispersioncorrect)
    obj.get_full_energy(rvdw=rvdw, rcoulomb=rcoulomb, dispersioncorrect=dispersioncorrect)
    obj.calc_energycoh(total_mass, nmols)


    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)



# =============================================================================
if __name__ == "__main__":

    main_app()
import argparse
import utils
import os
import sys
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

    desc = """Calculate the polymer size (Rg, Ree, ...) from a MD trajectory.
    This is part of the polyanagro library"""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.",
                        action="store", required=True, default=None)

    group2 = parser.add_mutually_exclusive_group(required=True)
    group2.add_argument("--topo", dest="topo",
                        help="A topology file in tpr, data or pdb format.\n"
                             "tpr --> GROMACS, dat --> LAMMPS, pdb --> OTHERS",
                        action="store", metavar="TPR|DATA|PDB")

    parser.add_argument("--stride", dest="stride", type=int,
                        help="Take a frame each stride frames, for example 10",
                        action="store", required=False, default=1)

    parser.add_argument("--fraction_trj_avg", dest="frac_avg", type=float,
                        help="""Fraction of the trajectory to calculate the averages. 
                        Example: 0.25 means that the 25%% first frames are discarted in the average calculation.""",
                        action="store", required=False, default=1.0)

    parser.add_argument("--e2e", dest="listee",
                        help="Calculate the end to end distances using the heads and tails of the chains."
                             "The format of the files must be a line for chain: ich ihead itail."
                             "The index must start in zero.",
                        action="store", metavar="FILE_WITH_DATA_End2EndAtoms", required=False, default=None)

    parser.add_argument("--e2acf", dest="e2acf",
                        help="Calculate the end to end autocorrelation function.",
                        action="store_true", required=False, default=False)

    parser.add_argument("--c2n", dest="listbb",
                        help="""Data to calculate the Cn of a polymer.
                             The parameter listbb can be either a pdb file template
                             (using the beta field 1:for backbone
                             and 0 for a branch atom) or a file with the following format:
                             i) a label [ mol01 ], ii) after a list of the backbone atoms in the mol01.
                             You need as much labels as chains or molecules in your system """)

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="pol_size.log")

    parser.add_argument("-d", "--distributions", dest="isdist",
                        help="Calculate Ree and Rg distributions",
                        action="store_true", required=False)

    parser.add_argument("--bondorientation", dest="isbondorientation",
                        help="Calculate intermolecular bond orientation",
                        action="store_true", required=False)

    parser.add_argument("--unwrap", dest="isunwrap", type=str2bool,
                        help="If 1 the coordinates provided are wrapped in the trajectory",
                        required=True, metavar="True or False, 1 or 0")

    parser.add_argument("--rg_massw", dest="isrgmass",
                        help="Calculate the mass weighted radius of gyration",
                        action="store_true", required=False)

    parser.add_argument("--isodf", dest="isodf",
                        help="Calculate 1st and 2nd Legendre polynomials for\n "
                             "the correlation between bonds in a polymer chain",
                        action="store_true", required=False)

    args = parser.parse_args()

    for itrj in args.traj:
        if not os.path.isfile(itrj):
            print("\nERROR: File {} does not exist\n".format(itrj))
            exit()

    if not os.path.isfile(args.topo):
        print("\nERROR: File {} does not exist\n".format(args.topo))
        exit()

    if args.listee and not os.path.isfile(args.listee):
        print("\nERROR: File {} does not exist\n".format(args.listee))
        exit()

    if args.listbb and not os.path.isfile(args.listbb):
        print("\nERROR: File {} does not exist\n".format(args.listbb))
        exit()

    if args.frac_avg > 1:
        print("\nERROR: Fraction must a number between 0 and 1\n".format(args.listbb))
        exit()

    return args

# =============================================================================
def get_listend(args):

    if args.listee is None:
        isree = False
        isreeacf = False
        listend2end = None
    else:
        isree = True
        if args.e2acf:
            isreeacf = True
        else:
            isreeacf = False
        listend2end = []
        with open(args.listee, 'r') as f:
            lines = f.readlines()
            for iline in lines:
                if iline.find("#") != -1:
                    continue
                # check line
                try:
                    ich = int(iline.split()[0])
                    ihead = int(iline.split()[1])
                    itail = int(iline.split()[2])
                    listend2end.append([ich, ihead, itail])
                except IndexError:
                    pass

    return isree, isreeacf, listend2end

# =============================================================================
def get_backbone_atoms(args, natoms, iatch):

    if args.listbb is None:
        iscn = False
        backbone_list_atoms = None
        isbbatom = None
    else:
        backbone_list_atoms = []
        isbbatom = natoms*[False]
        iscn = True
        fnamepath = args.listbb
        ext = os.path.splitext(fnamepath)[-1]
        ich = -1
        # dat files from REPLICATE_POLYMER starts at 0
        if ext == ".dat":
            with open(fnamepath, "r") as f:
                lines = f.readlines()
                for iline in lines:
                    if iline.find("[") != -1:
                        ich += 1
                        backbone_list_atoms.append([])
                    else:
                        if len(iline.replace(" ", "")) != 0:
                            n = iline.split()
                            for item in n:
                                backbone_list_atoms[ich].append(int(item))
                                isbbatom[int(item)] = True
        # ndx files from GROMACS start at 1
        elif ext == ".ndx":
            with open(fnamepath, "r") as f:
                lines = f.readlines()
                for iline in lines:
                    if iline.find("[") != -1:
                        ich += 1
                        backbone_list_atoms.append([])
                    else:
                        ll = iline.split()
                        for item in ll:
                            backbone_list_atoms[ich].append(int(item)-1)
                            isbbatom[int(item)-1] = True
        elif ext == ".pdb":
            with open(fnamepath, "r") as f:
                lines = f.readlines()
                for iline in lines:
                    if iline.find("ATOM") != -1 or iline.find("HETATM") != -1 :
                        iatom = int(iline[6:11])
                        occupancy = float(iline[54:60])
                        tempfactor = float(iline[60:66])
                        if tempfactor > 0.0:
                            isbbatom[int(iatom)-1] = False
                        else:
                            ich = iatch[int(iatom)-1]
                            isbbatom[int(iatom)-1] = True
                            try:
                                backbone_list_atoms[ich].append(int(iatom) - 1)
                            except IndexError:
                                backbone_list_atoms.append([])
                                backbone_list_atoms[ich].append(int(iatom) - 1)
        else:
            print("ERROR: ext {} is unknown in polymer_size.py (get_backbone_atoms())".format(ext))
            exit()

    return iscn, backbone_list_atoms, isbbatom

# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
                         Polymer size calculations 
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
    m += "\t\t\tpolymer_size".format(os.path.split(sys.argv[0])[1])
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
    # Create object to calculate
    objcalc = pag.Chain_Statistics(trj, dt=trj.dt, stride=args.stride, log=log)
    # Check end2end and backbone lists
    isree, isreeacf, listend2end = get_listend(args)
    iscnn, backbone_list_atoms, isbbatom = get_backbone_atoms(args, trj.topology.natoms, trj.topology._iatch)
    # Uwrap or not coordinates
    isunwrap = args.isunwrap
    # Calculate chain dimensions
    objcalc.calculate(listend2end, diroutput="./", isree=isree, isrg=True, iscn=iscnn, acfE2E=isreeacf,
                      distributions=args.isdist, molecularweight=True, calc_Cn_bonds_distances=True,
                      single_Cn_unitvector=False, begin=0, unwrap_pbc=isunwrap,
                      backbone_list_atoms=backbone_list_atoms, isbondorientation=args.isbondorientation,
                      isbbatom=isbbatom, isinternalchaindist=iscnn, isrgmass=args.isrgmass, isodf=args.isodf)

    objcalc.statistics(args.frac_avg)

    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)


# =============================================================================
if __name__ == "__main__":

    main_app()

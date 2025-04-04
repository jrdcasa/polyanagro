import argparse
import utils
import os
import sys
import datetime
import topology
import MDAnalysis as mda
from collections import defaultdict

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
    subparser = parser.add_subparsers(dest='command', required=True)

    gen = subparser.add_parser('generate')
    calc = subparser.add_parser('calculate')

    gen.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.",
                        action="store", required=True, default=None)

    gen.add_argument("--topo", dest="topo",
                        help="A topology file in tpr, data or pdb format.\n"
                             "tpr --> GROMACS, dat --> LAMMPS, pdb --> OTHERS",
                        action="store", metavar="TPR|DATA|PDB")

    gen.add_argument("--listbb", dest="listbb",
                        help="""Data to assign backbone atoms of a polymer.\n
                             The parameter listbb can be:\n
                             a) a pdb file template (using the beta field 1:for backbone
                             and 0 for a branch atom)\n 
                             b) a file with the following format:\n
                             i) a label [ mol01 ], ii) after a list of the backbone atoms in the mol01.\n
                             You need as much labels as chains or molecules in your system.\n
                             Keep in mind that in this case the indices start at 0. \n
                             c) all (--listbb all) to label all atoms as backbone ones.""",
                     action="store", required=True, default=None)

    gen.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="pol_bonddist.log")

    calc.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.",
                        action="store", required=True, default=None)

    calc.add_argument("--topo", dest="topo",
                        help="A topology file in tpr, data or pdb format.\n"
                             "tpr --> GROMACS, dat --> LAMMPS, pdb --> OTHERS",
                        action="store", metavar="TPR|DATA|PDB", required=True)

    calc.add_argument("-b", "--bonddist", dest="bonddistlist", nargs="+",
                      help="A list of labels contained in the file bonds_data_dist.ndx",
                      action="store", required=False, metavar="C2-C2_bb")

    calc.add_argument("-a", "--angdist", dest="angdistlist", nargs="+",
                        help="A list of labels contained in the file angles_data_dist.ndx",
                        action="store", required=False, metavar="C2-C2-C2")

    calc.add_argument("-d", "--dihdist", dest="dihdistlist", nargs="+",
                        help="A list of labels contained in the file dihedrals_data_dist.ndx",
                        action="store", required=False, metavar="C2-C2-C2-C2")

    calc.add_argument("-i", "--imprdist", dest="impdistlist", nargs="+",
                        help="A list of labels contained in the file impropers_data_dist.ndx",
                        action="store", required=False, metavar="C1-C2-C2-C2")

    calc.add_argument("--stride", dest="stride",
                        help="Take a frame each stride frames, for example 10",
                        action="store", required=False, default=1)

    calc.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="pol_bonddist_gen.log")

    calc.add_argument("--unwrap", dest="isunwrap", type=str2bool,
                        help="If True the coordinates provided are unwrapped",
                        required=True, metavar="True or False, 1 or 0")

    args = parser.parse_args()

    for itrj in args.traj:
        if not os.path.isfile(itrj):
            print("\nERROR: File {} does not exist\n".format(itrj))
            exit()

    if not os.path.isfile(args.topo):
        print("\nERROR: File {} does not exist\n".format(args.topo))
        exit()

    # if args.listbb and not os.path.isfile(args.listbb):
    #     print("\nERROR: File {} does not exist\n".format(args.listbb))
    #     exit()

    # if args.frac_avg > 1:
    #     print("\nERROR: Fraction must a number between 0 and 1\n".format(args.listbb))
    #     exit()

    return args

# =============================================================================
def get_backbone_atoms(args, nmols, natoms, iatch, logger=None):

    if args.listbb is None:
        backbone_list_atoms = None
        isbbatom = None
    elif args.listbb.upper() == "ALL":
        backbone_dict_atoms = defaultdict(list)
        backbone_list_atoms = []
        isbbatom = natoms * [False]
        for iatom, ich in enumerate(iatch):
            backbone_dict_atoms[ich].append(iatom)
            isbbatom[iatom] = "True"
        for ich, ivalues in backbone_dict_atoms.items():
            backbone_list_atoms.append(ivalues)
    else:
        backbone_list_atoms = []
        isbbatom = natoms*[False]
        fnamepath = args.listbb
        ext = os.path.splitext(fnamepath)
        ich = -1
        if ext[1] == ".dat":
            try:
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
                                    isbbatom[int(item)] = "True"
            except FileNotFoundError:
                m = "ERROR!!: File {} does not exist. Aborting ...".format(fnamepath)
                print(m) if logger is None else logger.info(m)
                exit()
        elif ext[1] == ".ndx":
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
                            isbbatom[int(item)-1] = "True"
        else:   # PDB File
            for imol in range(0, nmols):
                backbone_list_atoms.append([])
            isbbatom = natoms * [False]
            fnamepath = args.listbb
            try:
                with open(fnamepath, "r") as f:
                    lines = f.readlines()
                    for iline in lines:
                        if iline.find("ATOM") != -1 or iline.find("HETATM") != -1 :
                            item = int(float(iline[62:67]))   # Beta field in the PDB
                            iat = int(iline[6:11])-1
                            ich = iatch[iat]
                            if item == 0:
                                backbone_list_atoms[ich].append(iat)
                                isbbatom[iat] = "True"
            except FileNotFoundError:
                m = "ERROR!!: File {} does not exist. Aborting ...".format(fnamepath)
                print(m) if logger is None else logger.info(m)
                exit()

    return backbone_list_atoms, isbbatom

# =============================================================================
def print_header(version, logger_log=None):

    msg = """
    ***********************************************************************
            Bonded Distributions (bond, angle, dihedral, improper)
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
    m += "\t\t\tbonded_distribution".format(os.path.split(sys.argv[0])[1])
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

    # # Check backbone lists
    # nmols = len(trj.topology._nmols)
    # natoms = trj.topology.natoms
    # iatch = trj.topology._iatch
    # backbone_list_atoms, isbbatom = get_backbone_atoms(args, nmols, natoms, iatch, logger=log)
    # trj.topology._isbackbone = isbbatom

    if args.command == "generate":
        # Check backbone lists
        nmols = len(trj.topology._nmols)
        natoms = trj.topology.natoms
        iatch = trj.topology._iatch
        backbone_list_atoms, isbbatom = get_backbone_atoms(args, nmols, natoms, iatch, logger=log)
        trj.topology._isbackbone = isbbatom
        # Create object to generate bond, angles, dihedrals and impropers data
        objdist = pag.BondedDistributions(trj, dt=trj.dt, log=log)
        objdist.generate()
    elif args.command == "calculate":
        # Uwrap or not coordinates
        isunwrap = args.isunwrap
        bonddistlist = args.bonddistlist
        # Bond distribution =========
        if bonddistlist is not None:
            for item in bonddistlist:
                ndxfilename = "bonds_data_dist.ndx"
                objdist = pag.BondedDistributions(trj, dt=trj.dt, log=log)
                objdist.calculate(begin=0, typelabel="bond", unwrap_pbc=isunwrap,
                                  ndx_filename = ndxfilename,
                                  dist_name=item)
                del objdist
        # Angle distribution =========
        angdistlist = args.angdistlist
        if angdistlist is not None:
            for item in angdistlist:
                ndxfilename = "angle_data_dist.ndx"
                objdist = pag.BondedDistributions(trj, dt=trj.dt, log=log)
                objdist.calculate(begin=0, typelabel="angle", unwrap_pbc=isunwrap,
                                  ndx_filename = ndxfilename,
                                  dist_name=item)
                del objdist
        # Dihedral distribution =========
        dihdistlist = args.dihdistlist
        dihdistneigh = True
        if dihdistlist is not None:
            for item in dihdistlist:
                ndxfilename = "dihedral_data_dist.ndx"
                objdist = pag.BondedDistributions(trj, dt=trj.dt, log=log)
                objdist.calculate(begin=0, typelabel="dihedral", unwrap_pbc=isunwrap,
                                  ndx_filename = ndxfilename,
                                  dist_name=item, dihdistneigh=dihdistneigh)

                del objdist

        # Improper distribution ==========
        impdistlist = args.impdistlist
        if impdistlist is not None:
            for item in impdistlist:
                ndxfilename = "improper_data_dist.ndx"
                objdist = pag.BondedDistributions(trj, dt=trj.dt, log=log)
                objdist.calculate(begin=0, typelabel="improper", unwrap_pbc=isunwrap,
                                  ndx_filename = ndxfilename,
                                  dist_name=item)
                del objdist


    else:
        exit()

    # # Uwrap or not coordinates
    # isunwrap = args.isunwrap
    # # Bond distribution
    # objcalc.calculate(unwrap_pbc=isunwrap)
#     objcalc.calculate(listend2end, diroutput="./", isree=isree, isrg=True, iscn=iscnn, acfE2E=isreeacf,
#                       distributions=args.isdist, molecularweight=True, calc_Cn_bonds_distances=True,
#                       single_Cn_unitvector=False, begin=0, unwrap_pbc=isunwrap,
#                       backbone_list_atoms=backbone_list_atoms, isbondorientation=args.isbondorientation,
#                       isbbatom=isbbatom, isinternalchaindist=iscnn, isrgmass=args.isrgmass, isodf=args.isodf)
#
#     objcalc.statistics(args.frac_avg)
#
    now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    m = "\t\tJob  Done at {} ============\n".format(now)
    print(m) if log is None else log.info(m)
#
#
# =============================================================================
if __name__ == "__main__":

    main_app()

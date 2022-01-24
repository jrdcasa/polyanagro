import polyanagro as pag
import argparse
import topology
import utils
import os

# =============================================================================
def parse_arguments():

    desc = """Calculate the polymer size (Rg, Ree, ...) from a MD trajectory.
    This is part of the polyanagro library"""

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-t", "--traj", dest="traj", nargs="+",
                        help="A list of trajectories from MD simulations.",
                        action="store", required=True, default=None)

    parser.add_argument("--stride", dest="stride",
                        help="Take a frame each stride frames, for example 10",
                        action="store", required=False, default=1)

    parser.add_argument("--e2e", dest="listee",
                        help="Calculate the end to end distances using the heads and tails of the chains."
                             "The format of the files must be a line for chain: ich ihead itail."
                             "The index must start in zero.",
                        action="store", metavar="FILE_WITH_DATA", required=False, default=None)

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

    group2 = parser.add_mutually_exclusive_group(required=True)

    group2.add_argument("--tpr", dest="topo",
                        help="A topology file in tpr format.",
                        action="store")
    group2.add_argument("--psf", dest="topo",
                        help="A topology file in psf format.",
                        action="store")

    parser.add_argument("--log", dest="log",
                        help="Name of the file to write logs from this command",
                        action="store", required=False, default="pol_size.log")

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
def get_backbone_atoms(args, natoms):


    if args.listbb is None:
        iscn = False
        backbone_list_atoms = None
        isbbatom = None
    else:
        backbone_list_atoms = []
        isbbatom = natoms*[False]
        iscn = True
        fnamepath = args.listbb
        ext = os.path.splitext(fnamepath)
        ich = -1
        if ext != ".pdb":
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



    return iscn, backbone_list_atoms, isbbatom


# =============================================================================
def main_app():

    # Parse arguments
    args = parse_arguments()
    # Setup log
    log = utils.init_logger("Output", fileoutput=args.log, append=False, inscreen=False)
    # Load trajectory
    trj = topology.ExtTrajectory(args.traj, topfile=args.topo, logger=log)
    # Create object to calculate
    objcalc = pag.Chain_Statistics(trj, dt=trj.dt, stride=args.stride, log=log)
    # Check end2end and the list
    backbone_list_atoms = None
    is_bb_atoms = None

    isree, isreeacf, listend2end = get_listend(args)
    iscnn, backbone_list_atoms, isbbatom = get_backbone_atoms(args, trj.topology.natoms)

    objcalc.calculate(listend2end, diroutput="./", isree=isree, isrg=True, iscn=True, acfE2E=isreeacf,
                      distributions=False, molecularweight=True, calc_Cn_bonds_distances=True,
                      single_Cn_unitvector=True, begin=0, unwrap_pbc=True,
                      backbone_list_atoms=backbone_list_atoms,
                      isbbatom=isbbatom)




# =============================================================================
if __name__ == "__main__":

    main_app()
import polyanagro as pag
import numpy as np

# ===============================================================
class RDF(pag.Calculations):

    # #########################################################################
    def __init__(self, trj, iniframe=0, endframe=None, dt=1, stride=1,
                 setA_label="Full", setB_label="Full", onlybackbone=True, logger=None):

        super().__init__(trj, dt=dt, stride=stride, logger=logger)

        # Log file
        if not logger is None:
            self.logger = logger
        else:
            self.logger = None
        # Flag to know if the RDF is OK
        self._isok = True

        # Sets
        self._setA = None
        self._setB = None
        self._onlybackbone = onlybackbone

        # Define the sets to calculate the RDF
        self._isok = self._defineSets(setA_label, setB_label)

        # Frames
        self._iniframe = iniframe
        if endframe is None:
            self._endframe = self._trajectory.get_numframes()
        else:
            self._endframe = endframe
        self._stride = stride

        # Get the maximum distance in the the box dimensions
        self._max_box = -1.0
        for iframe in range(self._iniframe, self._endframe, self._stride):
            box_dimensions = self._trajectory.universe.trajectory[iframe].dimensions[0:3]
            m = max(box_dimensions)
            if m > self._max_box:
                self._max_box = m

    # #########################################################################
    def rdf_calc_cython(self):

        # Initialize the RDF arrays
        pag.setup_rdf_init(self._max_box)

        nat_A = self._setA.n_atoms
        nat_B = self._setB.n_atoms
        # atindex_A = np.array(self._setA.indices, dtype=np.int32)
        # b = np.array(self._setB.indices, dtype=np.int32)


        # Accumulate histograms
        for iframe in range(self._iniframe, self._endframe, self._stride):
            coords_t0_wrapped = self._trajectory.universe.trajectory[iframe].positions
            box_dimensions = self._trajectory.universe.trajectory[iframe].dimensions[0:3]


            isok = pag.rdf_calc(nat_A, nat_B, self._setA.indices, self._setB.indices,
                                self._setA.positions, self._setB.positions,
                                box_dimensions)

            print("isOK: {}", isok)

    # #########################################################################
    def _defineSets(self, setA_label, setB_label):

        # Take some info from the trajectory
        nmols_array, neigbors = self._trajectory.topology.get_array_mols_neigh()
        nchains = len(nmols_array)
        all_bb_bonds, nbonds_bb_perch = self._trajectory.topology.get_all_bb_bonds()

        # Check backbones calculations
        if self._onlybackbone:
            if len(self._trajectory.topology._isbackbone) != 0:
                str_atoms = ""
                for iat in range(self._trajectory.topology.natoms):
                    if self._trajectory.topology._isbackbone[iat]:
                        str_atoms += "{0:d} ".format(iat)
            else:
                str_atoms = None
                msg  = "\tRDF calculation for backbone cannot be performed\n"
                msg += "\tBackbone atoms have not been defined. \n"
                msg += "\tTry to use psf file as topology file. \n"
                print(msg) if self.logger is None else self.logger.warning(msg)
                return False
        else:
            str_atoms = ""
            for iat in range(self._trajectory.topology.natoms):
                str_atoms += "{0:d} ".format(iat)

        # Defining set groups
        if setA_label.upper() == "FULL" and setB_label.upper() == "FULL":
            self._setA = self._trajectory.universe.select_atoms("index {0:s}\n".format(str_atoms))
            self._setB = self._setA
        elif setA_label.upper() == "FULL" and not setB_label.upper() == "FULL":
            self._setA = self._trajectory.universe.select_atoms("index {0:s}\n".format(str_atoms))
            self._setB = self._trajectory.universe.select_atoms("{0:s}\n".format(setB_label))
        else:
            self._setA = self._trajectory.universe.select_atoms("{0:s}\n".format(setA_label))
            self._setB = self._trajectory.universe.select_atoms("{0:s}\n".format(setB_label))

        return True












        # Select both groups for the chain by index=========
        # for ich in range(nchains):
        #     str_atoms = ""
        #     for iat in nmols_array[ich]:
        #         str_atoms += "{0:d} ".format(iat)
        #     g1 = self._trajectory.universe.select_atoms("index {0:s}\n".format(str_atoms))
        #     g2 = self._trajectory.universe.select_atoms("index {0:s}\n".format(str_atoms))
        #     g1 = self._trajectory.universe.select_atoms("index {0:s}\n".format("0"))
        #     g2 = self._trajectory.universe.select_atoms("index {0:s}\n".format("2"))
        #     set1 = np.array(g1.positions, dtype=np.float64)
        #     set2 = np.array(g2.positions, dtype=np.float64)
        #     a = pag.distance_array(set1, set2)
        #
        #     print(a)
        #     print("=====================")


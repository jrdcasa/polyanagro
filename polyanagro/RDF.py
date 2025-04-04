import MDAnalysis.exceptions
import datetime
import polyanagro as pag
import numpy as np
import math


# ===============================================================
class RDF(pag.Calculations):

    #########################################################################
    def __init__(self, trj, iniframe=0, endframe=None, dt=1, stride=1, excl=None,
                 setA=None, setB=None, onlybackbone=True, delta_r = 0.02,
                 cutoff=None, logger=None):

        """

        Args:
            trj:
            iniframe:
            endframe:
            dt:
            stride:
            setA:
            setB:
            onlybackbone:
            delta_r: Delta r for the histogram in ansgtroms
            logger:
        """

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
        self._avg_volume = 0.0

        # Define the sets to calculate the RDF
        #self._isok = self._defineSets(setA_label, setB_label)

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

        self._delta_r = delta_r

        # Number of bins
        self._cutoff = cutoff
        if self._cutoff is None:
            maxDist = self._max_box * 0.5   # Half Box
        else:
            maxDist = self._cutoff
        self._nbins = int(maxDist / self._delta_r)

        m = "\t\t*********   RDF INFO   *********\n"
        m += "\t\t  Number of bins = {0:d}\n".format(self._nbins)
        m += "\t\t  Cutoff for RDF = {0:.2f} angstroms\n".format(self._cutoff)
        m += "\t\t  Bin width      = {0:.2f} angstroms\n".format(self._delta_r)
        m += "\t\t  Intramolecular = {0:b}\n".format(True)
        if excl is None:
            m += "\t\t  Exclusions     = None\n"
        else:
            m += "\t\t  Exclusions     = {0:d}\n".format(excl)
        m += "\t\t********* END RDF INFO *********\n"

        # Initialize the RDF arrays
        self._hist_total_rdf = np.zeros(self._nbins, dtype=np.int32)
        self._valid_counts_full_rdf = np.zeros(self._nbins, dtype=np.int32)
        self._hist_intra_rdf = np.zeros(self._nbins, dtype=np.int32)
        self._hist_inter_rdf = np.zeros(self._nbins, dtype=np.int32)
        self._hist_self_rdf = np.zeros(self._nbins, dtype=np.int32)

        self._total_rdf = np.zeros(self._nbins, dtype=np.float32)
        self._npairs_rdf_total = 0

        self._intra_rdf = np.zeros(self._nbins, dtype=np.float32)
        self._npairs_rdf_intra = 0

        self._inter_rdf = np.zeros(self._nbins, dtype=np.float32)
        self._npairs_rdf_inter = 0
        
        print(m) if self._logger is None else logger.info(m)

        # Exclude bonded interactions if 1 (1,2), if 2 (1,3), if 3 (1,4)
        self._maxnumberexclperatom = 40
        self._excl = excl
        # self._excl_array = np.zeros([self._trajectory.topology.natoms,
        #                              self._maxnumberexclperatom], dtype=np.int32)
        self._excl_array = np.full([self._trajectory.topology.natoms,
                                     self._maxnumberexclperatom], -1,  dtype=np.int32)
        if self._excl is not None:
            self._setup_exclusions()

    # #########################################################################
    def rdf_calc_cython(self):

        # Setup sets
        try:
            nat_A = len(self._setA)
            nat_B = len(self._setB)
        except TypeError:
            m = "\n\t\t Error!!! Set A and/or Set B are empty!!!\n"
            m += "\t\t          Set A : {}\n".format(self._setA)
            m += "\t\t          Set B : {}\n".format(self._setB)
            print(m) if self._logger is None else self._logger.error(m)
            exit()

        setA = np.array(self._setA, dtype=np.int32)
        setB = np.array(self._setB, dtype=np.int32)
        volumesum = 0
        nframes = 0

        # Accumulate histograms
        start_time = datetime.datetime.now()
        for iframe in range(self._iniframe, self._endframe, self._stride):

            if iframe % 100 == 0:
                msg = "\tRDF Frame {0:9d} of {1:9d}".format(iframe, self._endframe)
                mid_time = datetime.datetime.now()
                elapsed_time = mid_time - start_time
                msg += "\ttime: {0:s} seconds".format(str(elapsed_time.total_seconds()))
                print(msg) if self._logger is None else self._logger.info(msg)


            coords_t0_wrapped = self._trajectory.universe.trajectory[iframe].positions
            box_dimensions = self._trajectory.universe.trajectory[iframe].dimensions[0:3]
            volumesum += self._trajectory.universe.dimensions[0]*\
                         self._trajectory.universe.dimensions[1]*\
                         self._trajectory.universe.dimensions[2]

            if self._excl is None:
                self._npairs_rdf_total, \
                self._npairs_rdf_intra, \
                self._npairs_rdf_inter  = pag.rdf_hist_openmp(nat_A, nat_B,
                                                              self._nbins, self._delta_r, self._cutoff,
                                                              setA, setB,
                                                              coords_t0_wrapped,
                                                              box_dimensions,
                                                              self._trajectory.topology._iatch,
                                                              self._hist_total_rdf,
                                                              self._hist_intra_rdf,
                                                              self._hist_inter_rdf,
                                                              self._hist_self_rdf,
                                                             )
            else:
                self._npairs_rdf_total, \
                self._npairs_rdf_intra, \
                self._npairs_rdf_inter  = pag.rdf_hist_excl_openmp(nat_A, nat_B,
                                                                   self._nbins, self._delta_r, self._cutoff,
                                                                   setA, setB,
                                                                   self._excl_array,
                                                                   coords_t0_wrapped,
                                                                   box_dimensions,
                                                                   self._trajectory.topology._iatch,
                                                                   self._hist_total_rdf,
                                                                   self._hist_intra_rdf,
                                                                   self._hist_inter_rdf,
                                                                   self._hist_self_rdf,
                                                                  )
            nframes += 1

        self._avg_volume = volumesum / nframes
        # Get gr and normalizate
        isok = pag.rdf_gr(nframes, nat_A, nat_B, self._nbins,
                          self._cutoff, self._avg_volume,
                          self._hist_total_rdf, self._hist_intra_rdf,
                          self._hist_inter_rdf, self._hist_self_rdf,
                          self._total_rdf, self._npairs_rdf_total,
                          self._intra_rdf, self._npairs_rdf_intra,
                          self._inter_rdf, self._npairs_rdf_inter)

        self._rdf_write_to_file()


    # #########################################################################
    def define_sets(self, setA=None, setB=None):

        if setA is not None and setA.upper() == "NONE":
            setA = None
        if setB is not None and setB.upper() == "NONE":
            setB = None

        if setA is None and setB is None:
            self._setA = [i for i in range(0,self._trajectory.topology.natoms)]
            self._setB = [i for i in range(0, self._trajectory.topology.natoms)]
        elif setA is None and setB is not None:
            self._setA = [i for i in range(0,self._trajectory.topology.natoms)]
            self._setB = self._trajectory.universe.select_atoms("{0:s}\n".format(setB)).ids
        elif setA is not None and setB is None:
            self._setA = self._trajectory.universe.select_atoms("{0:s}\n".format(setA)).ids
            self._setB = [i for i in range(0,self._trajectory.topology.natoms)]
        elif setA is not None and setB is not None:
            try:
                self._setA = self._trajectory.universe.select_atoms("{0:s}\n".format(setA)).ids
            except MDAnalysis.exceptions.SelectionError:
                m = "\t\t\t setA label: {} \n".format(setA)
                m += "\t\t\t setA label is not allowed in MDAnalysis."
                print(m) if self._logger is None else self._logger.warn(m)
                self._setA = None
            try:
                self._setB = self._trajectory.universe.select_atoms("{0:s}\n".format(setB)).ids
            except MDAnalysis.exceptions.SelectionError:
                m = "\t\t\t setB label: {} \n".format(setB)
                m += "\t\t\t set B labels is not allowed in MDAnalysis."
                print(m) if self._logger is None else self._logger.warn(m)
                self._setB = None

        # # Take some info from the trajectory
        # nmols_array, neigbors = self._trajectory.topology.get_array_mols_neigh()
        # nchains = len(nmols_array)
        # all_bb_bonds, nbonds_bb_perch = self._trajectory.topology.get_all_bb_bonds()
        #
        # # Check backbones calculations
        # if self._onlybackbone:
        #     if len(self._trajectory.topology._isbackbone) != 0:
        #         str_atoms = ""
        #         for iat in range(self._trajectory.topology.natoms):
        #             if self._trajectory.topology._isbackbone[iat]:
        #                 str_atoms += "{0:d} ".format(iat)
        #     else:
        #         str_atoms = None
        #         msg  = "\tRDF calculation for backbone cannot be performed\n"
        #         msg += "\tBackbone atoms have not been defined. \n"
        #         msg += "\tTry to use psf file as topology file. \n"
        #         print(msg) if self.logger is None else self.logger.warning(msg)
        #         return False
        # else:
        #     str_atoms = ""
        #     for iat in range(self._trajectory.topology.natoms):
        #         str_atoms += "{0:d} ".format(iat)
        #
        # # Defining set groups
        # if setA_label.upper() == "FULL" and setB_label.upper() == "FULL":
        #     self._setA = self._trajectory.universe.select_atoms("index {0:s}\n".format(str_atoms))
        #     self._setB = self._setA
        # elif setA_label.upper() == "FULL" and not setB_label.upper() == "FULL":
        #     self._setA = self._trajectory.universe.select_atoms("index {0:s}\n".format(str_atoms))
        #     self._setB = self._trajectory.universe.select_atoms("{0:s}\n".format(setB_label))
        # else:
        #     self._setA = self._trajectory.universe.select_atoms("{0:s}\n".format(setA_label))
        #     self._setB = self._trajectory.universe.select_atoms("{0:s}\n".format(setB_label))

        self._writesetinfo()

        return True

    # #########################################################################
    def define_sets_idx(self, idx_file=None):

        setA = []
        setB = []
        iset = 0
        with open(idx_file, 'r') as fidx:
            lines = fidx.readlines()
            for iline in lines:
                if iline.find("[") == -1:
                    if iset == 1:
                        for token in iline.split():
                            setA.append(int(token)-1)
                    elif iset == 2:
                       for token in iline.split():
                            setB.append(int(token)-1)
                else:
                    iset += 1

        self._setA = setA
        self._setB = setB
        self._writesetinfo()


    # #########################################################################
    def _rdf_write_to_file(self, pattern=None):

        if pattern is None:
            filename_rdf_total = "rdf_total.dat"
            filename_rdf_intra = "rdf_intra.dat"
            filename_rdf_inter = "rdf_inter.dat"
        else:
            filename_rdf_total = "rdf_total_{}.dat".format(pattern)
            filename_rdf_intra = "rdf_intra_{}.dat".format(pattern)
            filename_rdf_inter = "rdf_inter_{}.dat".format(pattern)

        with open(filename_rdf_total, 'w') as frdf:

            for ibin in range(self._nbins):
                rl = ibin * self._delta_r
                ru = (ibin + 1) * self._delta_r
                rh = (rl + ru) / 2.0
                line = "{0:.2f} {1:.4f}\n".format(rh, self._total_rdf[ibin])
                frdf.writelines(line)

        with open(filename_rdf_intra, 'w') as frdf:

            for ibin in range(self._nbins):
                rl = ibin * self._delta_r
                ru = (ibin + 1) * self._delta_r
                rh = (rl + ru) / 2.0
                line = "{0:.2f} {1:.4f}\n".format(rh, self._intra_rdf[ibin])
                frdf.writelines(line)

        with open(filename_rdf_inter, 'w') as frdf:

            for ibin in range(self._nbins):
                rl = ibin * self._delta_r
                ru = (ibin + 1) * self._delta_r
                rh = (rl + ru) / 2.0
                line = "{0:.2f} {1:.4f}\n".format(rh, self._inter_rdf[ibin])
                frdf.writelines(line)

    # #########################################################################
    def _setup_exclusions(self):

        for iat in range(self._trajectory.topology.natoms):
            tmp_list = []
            for ilen in range(self._excl, 0, -1):
                paths = self._trajectory.topology.find_all_paths_length(iat, ilen)
                for ipath in paths:
                    for item in ipath:
                        if item == iat:
                            continue
                        tmp_list.append(item)
            tmp_set = set(tmp_list)
            if len(tmp_set) >= self._maxnumberexclperatom:
                m = "\n\t\tERROR!!!! Maximum number of neighbours exceeded\n"
                m += "\t\t\t  Maximum number of " \
                     "neighbours = {}\n".format(self._maxnumberexclperatom)
                m += "\t\t\t  iat = {}\n".format(iat)
                m += "\t\t\t  Increase the value of self._maxnumberexclpeatom in " \
                     "RDF.py or decrease the --excl number\n"
                print(m) if self._logger is None else self._logger.error(m)

            for icol, j in enumerate(tmp_set):
                self._excl_array[iat, icol] = j

    # #########################################################################
    def _writesetinfo(self):

        m = "\t\t*********   SET INFO   *********\n"
        m += "\t\t  Elements in setA = {0:d}\n".format(len(self._setA))
        m += "\t\t  Elements in setB = {0:d}\n".format(len(self._setB))
        if self._setA == self._setB:
            m += "\t\t  setA and setB are equals.\n".format(len(self._setB))
        else:
            m += "\t\t  setA and setB are differents.\n".format(len(self._setB))
        m += "\t\t********* END SET INFO *********\n"
        print(m) if self._logger is None else self._logger.info(m)

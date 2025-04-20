import datetime


import MDAnalysis as mda
import numpy as np
import re
import polyanagro as pag
from collections import defaultdict


class BondedDistributions(pag.Calculations):

    __slots__ = ["_deltaBond", "_maxbinBond", "_deltaAngle", "_maxbinAngle",
                 "_deltaDih", "_maxbinDih", "_deltaImp", "_maxbinImp", "_bdist",
                 "_adist", "_ddist", "_ddistFlory", "_idist", "_tacticitydist",
                 "_bondHist", "_angleHist", "_dihHist", "_dihHistFlory", "_impHist",
                 "_tactHist", "_filenameBondDist", "_fileAngleDist", "_fileDihDist",
                 "_fileImpDist", "_fileDihDistFlory", "_bond2DArray", "_angle2DList",
                 "_filenameAngleDist", "_dihedral2DArray", "_filenameDihedralDist",
                 "_improper2DArray", "_filenameImproperDist", "_dihvalues1DArray",
                 "_dihlabels1DArray", "_nmaxtList", "_nmaxgList", "_nmaxuList", "_keysUnits",
                 "_keysDyads", "_keysTryads", "_unitsDict", "_dyadsDict", "_tryadsDict", "_filenameDyadsDist",
                 "_dihedral2DExtendedArray", "_dihedral2DAllbb"]

    # #######################################################################
    def __init__(self, trj, dt=1, stride=1, log=None):

        super().__init__(trj, dt=dt, stride=stride, logger=log)

        #Histogram distribution characteristics
        self._deltaBond = 0.002
        self._maxbinBond = 2000

        self._deltaAngle  = 0.5
        self._maxbinAngle = 360

        self._deltaDih  = 1
        self._maxbinDih = 360

        self._deltaImp  = 1.0
        self._maxbinImp = 360

        #Distribution files
        self._filenameBondDist  = None
        self._filenameAngleDist  = None
        self._filenameDihedralDist = None
        self._filenameImproperDist = None

        #Bin Histograms
        self._bdist = np.zeros([self._maxbinBond],dtype=np.float64)
        self._adist = np.zeros([self._maxbinAngle],dtype=np.float64)
        self._ddist = np.zeros([self._maxbinDih],dtype=np.float64)
        self._ddistFlory = np.zeros([self._maxbinDih],dtype=np.float64)
        self._idist = np.zeros([self._maxbinImp],dtype=np.float64)
        self._tacticitydist = np.zeros([self._maxbinDih],dtype=np.float64)

        self._bondHist = np.zeros([self._maxbinBond],dtype=np.int32)
        self._angleHist = np.zeros([self._maxbinAngle],dtype=np.int32)
        self._dihHist = np.zeros([self._maxbinDih],dtype=np.int32)
        self._dihHistFlory = np.zeros([self._maxbinDih],dtype=np.int32)
        self._impHist = np.zeros([self._maxbinImp],dtype=np.int32)
        self._tactHist = np.zeros([self._maxbinDih],dtype=np.int32)

        self._bond2DArray = None
        self._angle2DList = None
        self._dihedral2DArray = None
        self._improper2DArray = None
        self._dihedral2DAllbb = None
        self._dihedral2DExtendedArray = None  # Take into account the next-neighbours

        self._dihvalues1DArray = None
        self._dihlabels1DArray = None

        #Setup histograms
        self._setupHistograms()

        self._nmaxtList = list()
        self._nmaxgList = list()
        self._nmaxuList = list()
        self._unitsDict = defaultdict(int)
        self._dyadsDict = defaultdict(int)
        self._tryadsDict = defaultdict(int)

        self._keysUnits = defaultdict()
        self._keysUnits['1'] = 'g'
        self._keysUnits['2'] = 't'
        self._keysUnits['3'] = 'u'

        self._keysDyads = defaultdict()
        self._keysDyads['11'] = 'gg'
        self._keysDyads['12'] = 'gt'
        self._keysDyads['13'] = 'gu'
        self._keysDyads['21'] = 'tg'
        self._keysDyads['22'] = 'tt'
        self._keysDyads['23'] = 'tu'
        self._keysDyads['31'] = 'ug'
        self._keysDyads['32'] = 'ut'
        self._keysDyads['33'] = 'uu'

        self._keysTryads = defaultdict()
        self._keysTryads['111'] = 'ggg'
        self._keysTryads['112'] = 'ggt'
        self._keysTryads['113'] = 'ggu'

        self._keysTryads['121'] = 'gtg'
        self._keysTryads['122'] = 'gtt'
        self._keysTryads['123'] = 'gtu'

        self._keysTryads['131'] = 'gug'
        self._keysTryads['132'] = 'gut'
        self._keysTryads['133'] = 'guu'

        self._keysTryads['211'] = 'tgg'
        self._keysTryads['212'] = 'tgt'
        self._keysTryads['213'] = 'tgu'

        self._keysTryads['221'] = 'ttg'
        self._keysTryads['222'] = 'ttt'
        self._keysTryads['223'] = 'ttu'

        self._keysTryads['231'] = 'tug'
        self._keysTryads['232'] = 'tut'
        self._keysTryads['233'] = 'tuu'

        self._keysTryads['311'] = 'ugg'
        self._keysTryads['312'] = 'ugt'
        self._keysTryads['313'] = 'ugu'

        self._keysTryads['321'] = 'utg'
        self._keysTryads['322'] = 'utt'
        self._keysTryads['323'] = 'utu'

        self._keysTryads['331'] = 'uug'
        self._keysTryads['332'] = 'uut'
        self._keysTryads['333'] = 'uuu'


    # #######################################################################
    def _setupHistograms(self):

        pag.setup_hist_bondC(self._deltaBond, self._maxbinBond, self._bdist)
        pag.setup_hist_angleC(self._deltaAngle, self._maxbinAngle, self._adist)
        pag.setup_hist_dihC(self._deltaDih, self._maxbinDih, self._ddist)
        pag.setup_hist_dihC(self._deltaDih, self._maxbinDih, self._ddistFlory)
        pag.setup_hist_dihC(self._deltaImp, self._maxbinImp, self._idist)
        pag.setup_hist_dihC(self._deltaDih, self._maxbinDih, self._tacticitydist)

    # #######################################################################
    def generate(self):

        now = datetime.datetime.now()
        m = "\n Generating the data to calculate bonded distributions ({})".format(now)
        print(m) if self._logger is None else self._logger.info(m)

        self._generate_bond_data()
        self._generate_angle_data()
        self._generate_dihedral_data()
        self._generate_improper_data()

        # Print labels in the log file
        try:
            with open("bonds_data_dist.ndx", 'r') as f:
                lines = f.readlines()
                pattern = "["
                matching_lines = [line for line in lines if pattern in line]
                m = "\t  * Bond types\n"
                m += "\t  ------------\n"
                for iline in matching_lines:
                    m += "\t  "+iline[1:-2]+"\n"
                m += "\t  ------------\n"
                print(m) if self._logger is None else self._logger.info(m)
        except FileNotFoundError:
            pass


        # Print labels in the log file
        try:
            with open("angle_data_dist.ndx", 'r') as f:
                lines = f.readlines()
                pattern = "["
                matching_lines = [line for line in lines if pattern in line]
                m = "\t  * Angles types\n"
                m += "\t  -------------\n"
                for iline in matching_lines:
                    m += "\t  "+iline[1:-2]+"\n"
                m += "\t  -------------\n"
                print(m) if self._logger is None else self._logger.info(m)
        except FileNotFoundError:
            pass

        # Print labels in the log file
        try:
            with open("dihedral_data_dist.ndx", 'r') as f:
                lines = f.readlines()
                pattern = "["
                matching_lines = [line for line in lines if pattern in line]
                m = "\t  * Dihedrals types\n"
                m += "\t  ----------------\n"
                for iline in matching_lines:
                    m += "\t  "+iline[1:-2]+"\n"
                m += "\t  ----------------\n"
                print(m) if self._logger is None else self._logger.info(m)
        except FileNotFoundError:
            pass

        # Print labels in the log file
        try:
            with open("improper_data_dist.ndx", 'r') as f:
                lines = f.readlines()
                pattern = "["
                matching_lines = [line for line in lines if pattern in line]
                m = "\t  * Improper types\n"
                m += "\t  ---------------\n"
                for iline in matching_lines:
                    m += "\t  "+iline[1:-2]+"\n"
                m += "\t  ---------------\n"
                print(m) if self._logger is None else self._logger.info(m)
        except FileNotFoundError:
            pass

    # #######################################################################
    def calculate(self, begin=0, unwrap_pbc=True, typelabel=None,
                  ndx_filename=None, dist_name=None, dihdistneigh=False):

        """
         Calculate bonded distributions

        Returns:

        """

        nframes = self._trajectory.get_numframes()
        if self._stride == 1:
            nframes_analysed = nframes
        else:
            nframes_analysed = int(((nframes - begin) / self._stride)) + 1
        m = "\t Num of frames to analyse: {}".format(nframes_analysed)
        print(m) if self._logger is None else self._logger.info(m)
        natoms = self._trajectory.natoms

        # Start calculations for each frame
        #nframes = self._trajectory.nframes
        ini = begin
        m = "\tProvide unwrapped PBC coordinates: {}".format(unwrap_pbc)
        print(m) if self._logger is None else self._logger.info(m)

        # Main loop of frames ====================================================
        s = datetime.datetime.now()
        idx_f = 0
        X1 = np.zeros([natoms], float)
        Y1 = np.zeros([natoms], float)
        Z1 = np.zeros([natoms], float)

        for iframe in range(ini, nframes, self._stride):

            # Estimated time (Use the 100 first frames to estimate the time)
            if idx_f == 0:
                f = datetime.datetime.now()
            if idx_f == 10:
                elapsed_time = datetime.datetime.now() - f
                k = int(((nframes - ini)/self._stride))/100
                estimated_time = k*elapsed_time.total_seconds()
                m = "\tESTIMATED TIME using {0:d} frames ({1:s} seconds): {2:.2f} seconds".format \
                    (100, str(elapsed_time.total_seconds()), estimated_time)
                print(m) if self._logger is None else self._logger.info(m)

            # Write info
            if iframe%self._freq == 0:
                elapsed_time = datetime.datetime.now() - s
                m = "\tIFRAME: {1:d} of {0:d} in {2:s} seconds".format \
                    (nframes_analysed, iframe, str(elapsed_time.total_seconds()))
                print(m) if self._logger is None else self._logger.info(m)

            # If pbc is false it is assumed that the trajectory is unwrapped
            if not unwrap_pbc:
                self._unwrap_coordinates(iframe)
            else:
                self._coords_unwrap = self._trajectory.universe.trajectory[iframe].positions

            X1[:] = self._coords_unwrap[:, 0]
            Y1[:] = self._coords_unwrap[:, 1]
            Z1[:] = self._coords_unwrap[:, 2]

            if typelabel == "bond":
                self._calculate_bond_dist_frame(ndx_filename, dist_name, X1, Y1, Z1)
            elif typelabel == "angle":
                self._calculate_angle_dist_frame(ndx_filename, dist_name, X1, Y1, Z1)
            elif typelabel == "dihedral":
                self._calculate_dihedral_dist_frame(ndx_filename, dist_name, X1, Y1, Z1, dihdistneigh=dihdistneigh)
            elif typelabel == "improper":
                self._calculate_improper_dist_frame(ndx_filename, dist_name, X1, Y1, Z1)
            elif typelabel == "dihneigh":
                self._calculate_dihneigh_dist_frame(ndx_filename, dist_name, X1, Y1, Z1)
            iframe += self._stride
            idx_f += 1

        self.writeDist(type=typelabel, ikey=dist_name)
        self._write_classify_torsions(dist_name)

        return True

    # #######################################################################
    def _generate_bond_data(self, issave=True):

        """
        Generate the necessary information to calculate the distribution of the bonds
        depending on the type of bond. The information is stored in several files.

        Returns:
            A list of all bonds, regardless of their type, and a dictionary listing
            bonds by type

        """

        # Get different type of bonds
        bond2DList = []
        typebond_dict = defaultdict(list)
        for ibond in self._trajectory.universe.bonds:
            at1 = ibond.indices[0]
            at2 = ibond.indices[1]
            itype1 = ibond.type[0]
            itype2 = ibond.type[1]

            # Equivalence labels
            if itype1 > itype2:
                strtype = itype1 + "-" + itype2
            else:
                strtype = itype2 + "-" + itype1
            typebond_dict[strtype].append([at1, at2])
            bond2DList.append([at1, at2])

        # Save a file for each type of bond
        bond_dict_tmp = defaultdict(list)
        if issave and len(typebond_dict) > 0:
            with open("bonds_data_dist.ndx", 'w') as fbond:
                fbond.writelines("# Index start at 1 (GROMACS). Bond types: {}\n".format(len(typebond_dict)))
                # Compile the bb and br atoms
                for ikey, ivalues in typebond_dict.items():
                    for ival in ivalues:
                        a = self._trajectory.topology._isbackbone[ival[0]]
                        b = self._trajectory.topology._isbackbone[ival[1]]
                        if a and b:
                            nkey = ikey + "_bb"
                        else:
                            nkey = ikey + "_br"
                        bond_dict_tmp[nkey].append(ival)

                # Write index bond file
                for ikey, ivalues in bond_dict_tmp.items():
                    fbond.writelines("[ {} ]\n".format(ikey))
                    for ival in ivalues:
                        fbond.writelines("{} {}\n".format(ival[0] + 1, ival[1] + 1))
                fbond.writelines("[ allbbbonds ]\n")
                for ikey, ivalues in bond_dict_tmp.items():
                    if ikey.find("bb") != -1:
                        for ival in ivalues:
                            fbond.writelines("{} {}\n".format(ival[0] + 1, ival[1] + 1))
                fbond.writelines("[ allbrbonds ]\n")
                for ikey, ivalues in bond_dict_tmp.items():
                    if ikey.find("br") != -1:
                        for ival in ivalues:
                            fbond.writelines("{} {}\n".format(ival[0] + 1, ival[1] + 1))
                fbond.writelines("[ allbonds ]\n")
                for ikey, ivalues in bond_dict_tmp.items():
                    for ival in ivalues:
                        fbond.writelines("{} {}\n".format(ival[0] + 1, ival[1] + 1))


        return bond2DList, typebond_dict

    # #######################################################################
    def _calculate_bond_dist_frame(self, ndx_filename, bonddist_label, X1, Y1, Z1):

        idx = -1
        if bonddist_label:
            with open(ndx_filename, 'r') as fbond:
                contents = fbond.readlines()
                fbond.seek(0)
                for num, line in enumerate(fbond, 1):
                    if re.match(r'.*'+bonddist_label+' ]$', line[:-1]):
                        idx = num
                # Only calculate the first frame
                if self._bond2DArray is None:
                    self._bond2DArray = []
                    if idx != -1:
                        while True:
                            try:
                                tmp = [int(i)-1 for i in contents[idx].split()]
                                self._bond2DArray.append(tmp)
                                idx += 1
                            except ValueError:
                                break
                            except IndexError:
                                m = "\t\t{} not label in {} index file".format(bonddist_label, ndx_filename)
                                print(m) if self._logger is None else self._logger.error(m)
                                exit()
                        self._bond2DArray = np.array(self._bond2DArray, dtype=np.int32)
                    else:
                        m = "\t\t{} not label in {} index file".format(bonddist_label, ndx_filename)
                        print(m) if self._logger is None else self._logger.error(m)
                        exit()

                iserror1 = self.bondDist(self._bond2DArray, X1, Y1, Z1)

    # #######################################################################
    def _generate_angle_data(self, issave=True):

        """
        Generate the necessary information to calculate the distribution of the angles
        depending on the type of angle. The information is stored in several files.

        Returns:
            A list of all angles, regardless of their type, and a dictionary listing
            angles by type

        """

        # # Get different type of bonds
        angle2DList = []
        typeangle_dict = defaultdict(list)
        for iangle in self._trajectory.universe.angles:
            at1 = iangle.indices[0]
            at2 = iangle.indices[1]
            at3 = iangle.indices[2]
            itype1 = iangle.type[0]
            itype2 = iangle.type[1]
            itype3 = iangle.type[2]

        # Equivalence labels
            if itype1 > itype3:
                strtype = itype1 + "-" + itype2 + "-" + itype3
            else:
                strtype = itype3 + "-" + itype2 + "-" + itype1
            typeangle_dict[strtype].append([at1, at2, at3])
            angle2DList.append([at1, at2, at3])

        # Save a file for each type of bond
        angle_dict_tmp = defaultdict(list)
        if issave and len(typeangle_dict) > 0:
            with open("angle_data_dist.ndx", 'w') as fangle:
                fangle.writelines("# Index start at 1 (GROMACS). Angle types: {}\n".format(len(typeangle_dict)))
                # Compile the bb and br atoms
                for ikey, ivalues in typeangle_dict.items():
                    for ival in ivalues:
                        a = self._trajectory.topology._isbackbone[ival[0]]
                        b = self._trajectory.topology._isbackbone[ival[1]]
                        c = self._trajectory.topology._isbackbone[ival[2]]
                        if a and b and c:
                            nkey = ikey + "_bb"
                        else:
                            nkey = ikey + "_br"
                        angle_dict_tmp[nkey].append(ival)

                # Write index angle file
                for ikey, ivalues in angle_dict_tmp.items():
                    fangle.writelines("[ {} ]\n".format(ikey))
                    for ival in ivalues:
                        fangle.writelines("{} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1))
                fangle.writelines("[ allbbangles ]\n")
                for ikey, ivalues in angle_dict_tmp.items():
                    if ikey.find("bb") != -1:
                        for ival in ivalues:
                            fangle.writelines("{} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1))
                fangle.writelines("[ allbrangles ]\n")
                for ikey, ivalues in angle_dict_tmp.items():
                    if ikey.find("br") != -1:
                        for ival in ivalues:
                            fangle.writelines("{} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1))
                fangle.writelines("[ allangles ]\n")
                for ikey, ivalues in angle_dict_tmp.items():
                    for ival in ivalues:
                        fangle.writelines("{} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1))

        self._angle2DList = angle2DList

        return angle2DList, typeangle_dict

    # #######################################################################
    def _calculate_angle_dist_frame(self, ndx_filename, angdist_label, X1, Y1, Z1):

        idx = -1
        if angdist_label:
            with open(ndx_filename, 'r') as fangle:
                contents = fangle.readlines()
                fangle.seek(0)
                for num, line in enumerate(fangle, 1):
                    if re.match(r'.*'+angdist_label+' ]$', line[:-1]):
                        idx = num
                # Only calculate the first frame
                if self._angle2DArray is None:
                    self._angle2DArray = []
                    if idx != -1:
                        while True:
                            try:
                                tmp = [int(i)-1 for i in contents[idx].split()]
                                self._angle2DArray.append(tmp)
                                idx += 1
                            except ValueError:
                                break
                            except IndexError:
                                m = "\t\t{} not label in {} index file".format(angdist_label, ndx_filename)
                                print(m) if self._logger is None else self._logger.error(m)
                                exit()
                        self._angle2DArray = np.array(self._angle2DArray, dtype=np.int32)
                    else:
                        m = "\t\t{} not label in {} index file".format(angdist_label, ndx_filename)
                        print(m) if self._logger is None else self._logger.error(m)
                        exit()

                iserror2 = self.angleDist(self._angle2DArray, X1, Y1, Z1)

    # #######################################################################
    def _generate_dihedral_data(self, issave=True):

        """
        Generate the necessary information to calculate the distribution of the angles
        depending on the type of angle. The information is stored in several files.

        Returns:
            A list of all angles, regardless of their type, and a dictionary listing
            angles by type

        """

        # Get different type of bonds
        dihedral2DList = []
        typedihedral_dict = defaultdict(list)
        for idhih in self._trajectory.universe.dihedrals:
            at1 = idhih.indices[0]
            at2 = idhih.indices[1]
            at3 = idhih.indices[2]
            at4 = idhih.indices[3]
            itype1 = idhih.type[0]
            itype2 = idhih.type[1]
            itype3 = idhih.type[2]
            itype4 = idhih.type[3]

            # Equivalence labels
            if itype1 == itype4:
                if itype2 > itype3:
                    strtype = itype1 + "-" + itype2 + "-" + itype3 + "-" + itype4
                else:
                    strtype = itype1 + "-" + itype3 + "-" + itype2 + "-" + itype4
            else:
                if itype1 > itype4:
                    if itype2 > itype3:
                        strtype = itype1 + "-" + itype2 + "-" + itype3 + "-" + itype4
                    else:
                        strtype = itype1 + "-" + itype3 + "-" + itype2 + "-" + itype4
                else:
                    if itype2 > itype3:
                        strtype = itype4 + "-" + itype2 + "-" + itype3 + "-" + itype1
                    else:
                        strtype = itype4 + "-" + itype3 + "-" + itype2 + "-" + itype1
            typedihedral_dict[strtype].append([at1, at2, at3, at4])
            dihedral2DList.append([at1, at2, at3, at4])

        # Save a file for each type of bond
        dihedral_dict_tmp = defaultdict(list)
        if issave and len(typedihedral_dict) > 0:
            with open("dihedral_data_dist.ndx", 'w') as fdih:
                fdih.writelines("# Index start at 1 (GROMACS). Dihedral types: {}\n".format(len(typedihedral_dict)))
                # Compile the bb and br atoms
                for ikey, ivalues in typedihedral_dict.items():
                    for ival in ivalues:
                        a = self._trajectory.topology._isbackbone[ival[0]]
                        b = self._trajectory.topology._isbackbone[ival[1]]
                        c = self._trajectory.topology._isbackbone[ival[2]]
                        d = self._trajectory.topology._isbackbone[ival[3]]
                        if a and b and c and d:
                            nkey = ikey + "_bb"
                        else:
                            nkey = ikey + "_br"
                        dihedral_dict_tmp[nkey].append(ival)

                # Write index angle file
                for ikey, ivalues in dihedral_dict_tmp.items():
                    fdih.writelines("[ {} ]\n".format(ikey))
                    for ival in ivalues:
                        fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))
                fdih.writelines("[ allbbangles ]\n")
                for ikey, ivalues in dihedral_dict_tmp.items():
                    if ikey.find("bb") != -1:
                        for ival in ivalues:
                            fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))
                fdih.writelines("[ allbrangles ]\n")
                for ikey, ivalues in dihedral_dict_tmp.items():
                    if ikey.find("br") != -1:
                        for ival in ivalues:
                            fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))
                fdih.writelines("[ allangles ]\n")
                for ikey, ivalues in dihedral_dict_tmp.items():
                    for ival in ivalues:
                        fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))

        return dihedral2DList, typedihedral_dict

    # #######################################################################
    def _calculate_dihedral_dist_frame(self, ndx_filename, dihdist_label, X1, Y1, Z1, dihdistneigh=False):

        idx = -1
        if dihdist_label:
            with open(ndx_filename, 'r') as fdih:
                contents = fdih.readlines()
                fdih.seek(0)
                for num, line in enumerate(fdih, 1):
                    if re.match(r'.*'+dihdist_label+' ]$', line):
                        idx = num
                # Only calculate the first frame
                if self._dihedral2DArray is None:
                    self._dihedral2DArray = []
                    if idx != -1:
                        while True:
                            try:
                                tmp = [int(i)-1 for i in contents[idx].split()]
                                self._dihedral2DArray.append(tmp)
                                idx += 1
                            except ValueError:
                                break
                            except IndexError:
                                m = "\t\t{} not label in {} index file".format(dihdist_label, ndx_filename)
                                print(m) if self._logger is None else self._logger.error(m)
                                exit()
                        self._dihedral2DArray = np.array(self._dihedral2DArray, dtype=np.int32)
                    else:
                        m = "\t\t{} not label in {} index file".format(dihdist_label, ndx_filename)
                        print(m) if self._logger is None else self._logger.error(m)
                        exit()

                iserror3 = self.dihDist(self._dihedral2DArray, X1, Y1, Z1, self._dihHist)

        # Get an aaray with all dihedral in the backbone sorted by the first column
        idx = -1
        if dihdistneigh:
            # Read allbbdihedrals:
            dihdist_label = "allbbangles"
            with open(ndx_filename, 'r') as fdih:
                fdih.seek(0)
                contents = fdih.readlines()
                fdih.seek(0)
                for num, line in enumerate(fdih, 1):
                    if re.match(r'.*' + dihdist_label + ' ]$', line):
                        idx = num

                # Only calculate the first frame
                if self._dihedral2DAllbb is None:
                    self._dihedral2DAllbb = []
                    if idx != -1:
                        while True:
                            try:
                                tmp = [int(i) - 1 for i in contents[idx].split()]
                                self._dihedral2DAllbb.append(tmp)
                                idx += 1
                            except ValueError:
                                break
                            except IndexError:
                                if idx == len(contents):
                                    break
                                else:
                                    m = "\t\t{} not label in {} index file".format(dihdist_label, ndx_filename)
                                    print(m) if self._logger is None else self._logger.error(m)
                                    exit()
                        self._dihedral2DAllbb = \
                            np.array(self._dihedral2DAllbb, dtype=np.int32)
                    else:
                        m = "\t\t{} not label in {} index file".format(dihdist_label, ndx_filename)
                        print(m) if self._logger is None else self._logger.error(m)
                        exit()

                    # Order the dihedrals using the first column
                    self._dihedral2DAllbb = self._dihedral2DAllbb[self._dihedral2DAllbb[:,0].argsort()]

            # # Loop over all dihedrals and look for the next neighbours
            # if self._dihedral2DExtendedArray is None:
            #     self._dihedral2DExtendedArray = []
            #     for item in self._dihedral2DArray:
            #
            #         idx_left_dihedral = np.where((self._dihedral2DAllbb[:, 1] == item[0]) &
            #                                      (self._dihedral2DAllbb[:, 2] == item[1]) &
            #                                      (self._dihedral2DAllbb[:, 3] == item[2]))[0]
            #         try:
            #             left_dihedral = self._dihedral2DAllbb[idx_left_dihedral, :][0]
            #             self._dihedral2DExtendedArray.append(left_dihedral)
            #
            #         except IndexError:
            #             left_dihedral = None
            #
            #         idx_right_dihedral = np.where((self._dihedral2DAllbb[:, 0] == item[1]) &
            #                                       (self._dihedral2DAllbb[:, 1] == item[2]) &
            #                                       (self._dihedral2DAllbb[:, 2] == item[3]))[0]
            #         try:
            #             right_dihedral = self._dihedral2DAllbb[idx_right_dihedral, :][0]
            #             self._dihedral2DExtendedArray.append(right_dihedral)
            #         except IndexError:
            #             right_dihedral = None
            #
            #         self._dihedral2DExtendedArray.append(item)
            #
            # self._dihedral2DExtendedArray = np.array(self._dihedral2DExtendedArray, dtype=np.int32)
            # self._dihvalues1DArray = np.zeros(len(self._dihedral2DExtendedArray), dtype=np.float64)
            # self._dihlabels1DArray = np.zeros(len(self._dihedral2DExtendedArray), dtype=np.int32)
            self._dihedral2DAllbb = np.array(self._dihedral2DAllbb, dtype=np.int32)
            self._dihvalues1DArray = np.zeros(len(self._dihedral2DAllbb), dtype=np.float64)
            self._dihlabels1DArray = np.zeros(len(self._dihedral2DAllbb), dtype=np.int32)
            iserror3 = self.dihDistNeigh(self._dihedral2DAllbb, X1, Y1, Z1,
                                         self._dihvalues1DArray, self._dihlabels1DArray)

            self._classify_torsions(self._dihedral2DAllbb, self._dihlabels1DArray)

    # #######################################################################
    def _generate_improper_data(self, issave=True):

        """
        Generate the necessary information to calculate the distribution of the angles
        depending on the type of angle. The information is stored in several files.

        Returns:
            A list of all angles, regardless of their type, and a dictionary listing
            angles by type

        """

        # Get different type of bonds
        improper2DList = []
        typeimproper_dict = defaultdict(list)
        for iimpr in self._trajectory.universe.impropers:
            at1 = iimpr.indices[0]
            at2 = iimpr.indices[1]
            at3 = iimpr.indices[2]
            at4 = iimpr.indices[3]
            itype1 = iimpr.type[0]
            itype2 = iimpr.type[1]
            itype3 = iimpr.type[2]
            itype4 = iimpr.type[3]

            strtype = itype1 + "-" + itype2 + "-" + itype3 + "-" + itype4
            typeimproper_dict[strtype].append([at1, at2, at3, at4])
            improper2DList.append([at1, at2, at3, at4])

        if len(typeimproper_dict) == 0:
            angles_bb = []
            for iangle in self._trajectory.universe.angles:
                at1 = iangle.indices[0]
                at2 = iangle.indices[1]
                at3 = iangle.indices[2]
                a = self._trajectory.topology._isbackbone[at1]
                b = self._trajectory.topology._isbackbone[at2]
                c = self._trajectory.topology._isbackbone[at3]
                if a and b and c :
                    angles_bb.append(iangle)
            ll = mda.topology.guessers.guess_improper_dihedrals(angles_bb)
            for iimpr in ll:
                at1 = iimpr[0]
                at2 = iimpr[1]
                at3 = iimpr[2]
                at4 = iimpr[3]
                itype1 = self._trajectory.universe.atoms[at1].type
                itype2 = self._trajectory.universe.atoms[at2].type
                itype3 = self._trajectory.universe.atoms[at3].type
                itype4 = self._trajectory.universe.atoms[at4].type

                strtype = itype1 + "-" + itype2 + "-" + itype3 + "-" + itype4
                typeimproper_dict[strtype].append([at1, at2, at3, at4])
                improper2DList.append([at1, at2, at3, at4])

        # Save a file for each type of bond
        improper_dict_tmp = defaultdict(list)
        if issave and len(typeimproper_dict) > 0:
            filename = "improper_data_dist.ndx"
            with open(filename, 'w') as fdih:
                fdih.writelines("# Index start at 1 (GROMACS). Improper types: {}\n".
                                format(len(typeimproper_dict)))
                # Compile the bb and br atoms
                for ikey, ivalues in typeimproper_dict.items():
                    for ival in ivalues:
                        a = self._trajectory.topology._isbackbone[ival[0]]
                        b = self._trajectory.topology._isbackbone[ival[1]]
                        c = self._trajectory.topology._isbackbone[ival[2]]
                        d = self._trajectory.topology._isbackbone[ival[3]]
                        if a and b and c and d:
                            nkey = ikey + "_bb"
                        else:
                            nkey = ikey + "_br"
                        improper_dict_tmp[nkey].append(ival)
                # Write index dihedral file
                for ikey, ivalues in typeimproper_dict.items():
                        fdih.writelines("[ {} ]\n".format(ikey))
                        for ival in ivalues:
                            fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))
                fdih.writelines("[ allbbimpropers ]\n")
                for ikey, ivalues in improper_dict_tmp.items():
                    if ikey.find("bb") != -1:
                        for ival in ivalues:
                            fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))
                fdih.writelines("[ allbrimpropers ]\n")
                for ikey, ivalues in improper_dict_tmp.items():
                    if ikey.find("br") != -1:
                        for ival in ivalues:
                            fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))
                fdih.writelines("[ allimpropers ]\n")
                for ikey, ivalues in improper_dict_tmp.items():
                    for ival in ivalues:
                        fdih.writelines("{} {} {} {}\n".format(ival[0] + 1, ival[1] + 1, ival[2] + 1, ival[3] + 1))

        return improper2DList, typeimproper_dict

    # #######################################################################
    def _calculate_improper_dist_frame(self, ndx_filename, imprdist_label, X1, Y1, Z1):

        idx = -1
        if imprdist_label:
            with open(ndx_filename, 'r') as fimpr:
                contents = fimpr.readlines()
                fimpr.seek(0)
                for num, line in enumerate(fimpr, 1):
                    if re.match(r'.*' + imprdist_label + ' ]$', line[:-1]):
                        idx = num
                # Only calculate the first frame
                if self._improper2DArray is None:
                    self._improper2DArray = []
                    if idx != -1:
                        while True:
                            try:
                                tmp = [int(i) - 1 for i in contents[idx].split()]
                                self._improper2DArray.append(tmp)
                                idx += 1
                            except ValueError:
                                break
                            except IndexError:
                                break
                                # m = "\t\t{} not label in {} index file".format(imprdist_label, ndx_filename)
                                # print(m) if self._logger is None else self._logger.error(m)
                                # exit()
                        self._improper2DArray = np.array(self._improper2DArray, dtype=np.int32)
                    else:
                        m = "\t\t{} not label in {} index file".format(imprdist_label, ndx_filename)
                        print(m) if self._logger is None else self._logger.error(m)
                        exit()

                iserror3 = self.imprDist(self._improper2DArray, X1, Y1, Z1)

    # #######################################################################
    def bondDist(self, bond2DArray, X1, Y1, Z1):

        iserror = 1
        if len(bond2DArray) != 0:
            iserror = pag.bondDistC(bond2DArray, X1, Y1, Z1, self._bondHist)
        return iserror

    ########################################################################
    def angleDist(self,angle2DArray,X1,Y1,Z1):

        iserror = 1
        if len(angle2DArray) != 0:
           iserror = pag.angleDistC(angle2DArray, X1, Y1, Z1, self._angleHist)
        return iserror


    ########################################################################
    @staticmethod
    def dihDist(dih2DArray, X1, Y1, Z1, dihHist):

        iserror = 1
        if len(dih2DArray) != 0:
            iserror = pag.dihDistC(dih2DArray, X1, Y1, Z1, dihHist)
        return iserror


    ########################################################################
    @staticmethod
    def dihDistNeigh(dih2DArray, X1, Y1, Z1, dihvalues1DArray, dihlabels1DAray):

        iserror = 1
        if len(dih2DArray) != 0:
            iserror = pag.dihDistCNeigh(dih2DArray, X1, Y1, Z1, dihvalues1DArray, dihlabels1DAray)
        return iserror

    ########################################################################
    def _classify_torsions(self, dihedralArray, dihlabelsArray):

        dictLabelTorsions = defaultdict(list)
        ind = 0

        for item in dihedralArray:

            ich1 = self._trajectory.topology._iatch[int(item[0])]
            ich2 = self._trajectory.topology._iatch[int(item[1])]
            ich3 = self._trajectory.topology._iatch[int(item[2])]
            ich4 = self._trajectory.topology._iatch[int(item[3])]
            l = [ich1, ich2, ich3, ich4]
            if l.count(l[0]) != len(l):
                m = "\n\t\t Dihedral {0} {1} {2}\n".format(item, l(l[0]), len(l))
                m += "\t\t Error. Atoms in dihedral not in the same chain"
                print(m) if self._logger is None else self._logger.error(m)
                exit()
            dictLabelTorsions[l[0]].append(str(dihlabelsArray[ind]))
            ind += 1

        # Dyads and triads
        for key in dictLabelTorsions.keys():
            dictLabelTorsions[key] = ''.join(dictLabelTorsions[key])
            try:
                self._nmaxtList.append(len(max(re.findall("2*2", dictLabelTorsions[key]))))
            except ValueError:
                self._nmaxtList.append(0)
            try:
                self._nmaxgList.append(len(max(re.findall("1*1", dictLabelTorsions[key]))))
            except ValueError:
                self._nmaxgList.append(0)
            try:
                self._nmaxuList.append(len(max(re.findall("3*3", dictLabelTorsions[key]))))
            except ValueError:
                self._nmaxuList.append(0)

            for i in range(0, len(dictLabelTorsions[key])):
                pair = dictLabelTorsions[key][i]
                label = self._keysUnits[pair]
                self._unitsDict[label] += 1

            for i in range(1, len(dictLabelTorsions[key])):
                pair = dictLabelTorsions[key][i - 1:i + 1]
                label = self._keysDyads[pair]
                self._dyadsDict[label] += 1

            for i in range(2, len(dictLabelTorsions[key])):
                pair = dictLabelTorsions[key][i - 2:i + 1]
                label = self._keysTryads[pair]
                self._tryadsDict[label] += 1

    #######################################################################
    def _write_classify_torsions(self, name):

        self._filenameDyadsDist = "Summary_Dyads_{}.dat".format(name)
        with open(self._filenameDyadsDist, "w") as fout:
            l1 = "<nmax>-trans  : {0:.2f} +- {1:.2f} max: {2:d} min: {3:d}\n". \
                format(np.mean(self._nmaxtList), np.std(self._nmaxtList), np.max(self._nmaxtList),
                       np.min(self._nmaxtList))
            l2 = "<nmax>-gauche : {0:.2f} +- {1:.2f} max: {2:d} min: {3:d}\n". \
                format(np.mean(self._nmaxgList), np.std(self._nmaxgList), np.max(self._nmaxgList),
                       np.min(self._nmaxgList))
            l3 = "<nmax>-gauche-: {0:.2f} +- {1:.2f} max: {2:d} min: {3:d}\n". \
                format(np.mean(self._nmaxuList), np.std(self._nmaxgList), np.max(self._nmaxuList),
                       np.min(self._nmaxuList))
            fout.writelines(l1)
            fout.writelines(l2)
            fout.writelines(l3)
            fout.writelines("\n\n")

            nt = 0
            for item in self._unitsDict:
                nt += self._unitsDict[item]

            fout.writelines("========= Units =============\n")
            for item in self._unitsDict:
                l = "{0:s}: {1:.1f} % \n".format(item, float(self._unitsDict[item]) / float(nt) * 100.)
                fout.writelines(l)

            fout.writelines("\n\n")

            nt = 0
            for item in self._dyadsDict:
                nt += self._dyadsDict[item]

            fout.writelines("========= Dyads =============\n")
            for item in self._dyadsDict:
                l = "{0:s}: {1:.1f} % \n".format(item, float(self._dyadsDict[item]) / float(nt) * 100.)
                fout.writelines(l)

            fout.writelines("\n\n")

            nt = 0
            for item in self._tryadsDict:
                nt += self._tryadsDict[item]
            fout.writelines("========= Tryads =============\n")
            for item in self._tryadsDict:
                l = "{0:s}: {1:.1f} % \n".format(item, float(self._tryadsDict[item]) / float(nt) * 100.)
                fout.writelines(l)


        totaldiah = sum([item for key, item in self._dyadsDict.items()])
        list_labels = ["tt", "gg", "uu", "tg", "tu"]
        line = ""
        idx = 0
        with open("Dist_Dyads_{}.dat".format(name), "w") as f:
            line = "0 tt {0:d} {1:7.5f}\n". \
                format(self._dyadsDict["tt"], float(self._dyadsDict["tt"] / float(totaldiah)))
            line += "1 gg {0:d} {1:7.5f}\n". \
                format(self._dyadsDict["gg"], float(self._dyadsDict["gg"] / float(totaldiah)))
            line += "2 uu {0:d} {1:7.5f}\n". \
                format(self._dyadsDict["uu"], float(self._dyadsDict["uu"] / float(totaldiah)))
            line += "3 tg {0:d} {1:7.5f}\n". \
                format(self._dyadsDict["tg"], float(self._dyadsDict["tg"] / float(totaldiah)))


            # line += ("3 tg %d %7.5f\n" % (self._dyadsDict('tg\n') + self._dyadsDict('gt\n'),
            #                              float(self._dyadsDict('tg\n') + self._dyadsDict('gt\n')) / float(totaldiah)))
            # line = ("4 tu %d %7.5f\n" % (self._dyadsDict('tu\n') + self._dyadsDict('ut\n'),
            #                              float(self._dyadsDict('tu\n') + self._dyadsDict('ut\n')) / float(totaldiah)))
            f.writelines(line)
            # line = ("5 gu %d %7.5f\n" % (self._dyadsDict('gu\n') + self._dyadsDict('ug\n'),
            #                              float(self._dyadsDict('gu\n') + self._dyadsDict('ug\n')) / float(totaldiah)))

    # ########################################################################
#     def dihDistFlory(self,dih2DArray,X1,Y1,Z1):
#
#         iserror = 1
#         if len(dih2DArray) != 0:
#             iserror = distC.dihDistFloryC(dih2DArray, X1, Y1, Z1, self._dihHistFlory)
#         return iserror
#
    # #######################################################################
    def imprDist(self,imp2DArray,X1,Y1,Z1):

        iserror = 1
        if len(imp2DArray) != 0:
            iserror = pag.dihDistC(imp2DArray, X1, Y1, Z1, self._impHist)
        return iserror

# ########################################################################
#     def tacticity(self,tact2DArray,X1,Y1,Z1):
#
#         iserror = 1
#         if len(tact2DArray) != 0:
#             iserror = distC.tactDistC(tact2DArray, X1, Y1, Z1, self._tactHist)
#         return iserror
#
# ########################################################################
#     @staticmethod
#     def write_tacticity_hist(self):
#
#         """
#                                    (at4)
#                                  D1  |    D3
#            (at7)---(at5)---(at2)---(at1) --- (at3)---(at6)
#
#            D1 (R1) --> 3-1-2-5
#            D2 (R2) --> 1-2-5-7
#            D3 (R3) --> 6-3-1-2
#         """
#
#         # Calculate the histogram for tacticity
#         f = open('tacticity_D.dat', 'r')
#         valdih1 = []
#         valdih2 = []
#         valdih3 = []
#         valdih1_3 = []
#         while 1:
#             line = f.readline()
#             if not line: break
#             d1 = line.split()[0]
#             d2 = line.split()[1]
#             d3 = line.split()[2]
#             valdih1.append(float(d1))
#             valdih2.append(float(d2))
#             valdih3.append(float(d3))
#             valdih1_3.append(float(d1))
#             valdih1_3.append(float(d3))
#         f.close()
#
#         # hist = np.histogram(valdih,bins=1,range=[-180,180],normed=True)
#         #freq, bins = np.histogram(valdih, bins=360, range=[-180, 180], normed=True)
#         freq1, bins1 = np.histogram(valdih1, bins=360, range=[0, 360], density=True)
#         data1 = zip(bins1, freq1)
#         freq2, bins2 = np.histogram(valdih2, bins=360, range=[0, 360], density=True)
#         data2 = zip(bins2, freq2)
#         freq3, bins3 = np.histogram(valdih3, bins=360, range=[0, 360], density=True)
#         data3 = zip(bins3, freq3)
#         freq1_3, bins1_3 = np.histogram(valdih1_3, bins=360, range=[0, 360], density=True)
#         data1_3 = zip(bins1_3, freq1_3)
#
#         f = open('Dih_Dist_tactD1.dat', 'w')
#         for item in data1:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#         f = open('Dih_Dist_tactD2.dat', 'w')
#         for item in data2:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactD3.dat', 'w')
#         for item in data3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactD1_3.dat', 'w')
#         for item in data1_3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         # Calculate the histogram for tacticity
#         f = open('tacticity_R.dat', 'r')
#         valdih1 = []
#         valdih2 = []
#         valdih3 = []
#         valdih1_3 = []
#         while 1:
#              line = f.readline()
#              if not line: break
#              d1 = line.split()[0]
#              d2 = line.split()[1]
#              d3 = line.split()[2]
#              valdih1.append(float(d1))
#              valdih2.append(float(d2))
#              valdih3.append(float(d3))
#              valdih1_3.append(float(d1))
#              valdih1_3.append(float(d3))
#
#         f.close()
#
#         # hist = np.histogram(valdih,bins=1,range=[-180,180],normed=True)
#         #freq, bins = np.histogram(valdih, bins=360, range=[-180, 180], normed=True)
#         freq1, bins1 = np.histogram(valdih1, bins=360, range=[0, 360], density=True)
#         data1 = zip(bins1, freq1)
#         freq2, bins2 = np.histogram(valdih2, bins=360, range=[0, 360], density=True)
#         data2 = zip(bins2, freq2)
#         freq3, bins3 = np.histogram(valdih3, bins=360, range=[0, 360], density=True)
#         data3 = zip(bins3, freq3)
#         freq1_3, bins1_3 = np.histogram(valdih1_3, bins=360, range=[0, 360], density=True)
#         data1_3 = zip(bins1_3, freq1_3)
#
#         f = open('Dih_Dist_tactR1.dat', 'w')
#         for item in data1:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactR2.dat', 'w')
#         for item in data2:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactR3.dat', 'w')
#         for item in data3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#         f = open('Dih_Dist_tactR1_3.dat', 'w')
#         for item in data1_3:
#             line= str(item[0])+" "+str(item[1])+"\n"
#             f.writelines(line)
#         f.close()
#
#
########################################################################
    def writeDist(self, type="all", ikey=None):

        if type == "bond":
            self._filenameBondDist = "Dist_Bond_{}.dat".format(ikey)
            with open(self._filenameBondDist, 'w') as f:
                for i in range(0,self._maxbinBond):
                    if np.sum(self._bondHist) == 0: break
                    f.writelines(str(self._bdist[i])+" "+str(self._bondHist[i])+" "+\
                                 str(self._bondHist[i]*(1.0/np.sum(self._bondHist))*(1.0/self._deltaBond))+"\n")
        if type == "angle":
            self._filenameAngleDist = "Dist_Angle_{}.dat".format(ikey)
            with open(self._filenameAngleDist, 'w') as f:
                for i in range(0,self._maxbinAngle):
                    if np.sum(self._angleHist) == 0: break
                    f.writelines(str(self._adist[i])+" "+str(self._angleHist[i])+" "+\
                                 str(self._angleHist[i]*(1.0/np.sum(self._angleHist))*(1.0/self._deltaAngle))+"\n")
        if type == "dihedral":
            self._filenameDihedralDist = "Dist_Dihedral_{}.dat".format(ikey)
            ddist_range = list()
            for i in self._ddist:
                if i >= 180:
                    ddist_range.append(int(i-360))
                else:
                    ddist_range.append(int(i))

            with open(self._filenameDihedralDist, 'w') as f:
                for i in range(0,self._maxbinDih):
                    if np.sum(self._dihHist) != 0:
                        line = "{0:>4d}  {1:>4d}  {2:>10d}  {3:>.3e}\n".\
                            format(int(self._ddist[i]), ddist_range[i], self._dihHist[i],
                                   self._dihHist[i]*(1.0/np.sum(self._dihHist))*(1.0/self._deltaDih) )
                        f.writelines(line)
        if type == "improper":
            self._filenameImproperDist = "Dist_Improper_{}.dat".format(ikey)
            with open(self._filenameImproperDist, 'w') as f:
                for i in range(0,self._maxbinImp):
                    if np.sum(self._impHist) != 0:
                        f.writelines(str(self._idist[i])+" "+str(self._impHist[i])+" "+\
                                         str(self._impHist[i]*(1.0/np.sum(self._impHist))*(1.0/self._deltaImp))+"\n")
        #
        # for i in range(0,self._maxbinDih):
        #     if np.sum(self._dihHistFlory) != 0:
        #         self._fileDihDistFlory.writelines(str(self._ddistFlory[i])+" "+str(self._dihHistFlory[i])+" "+\
        #                                           str(self._dihHistFlory[i]*(1.0/np.sum(self._dihHistFlory))*(1.0/self._deltaDih))+"\n")

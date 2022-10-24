import datetime
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
                 "_fileImpDist", "_fileDihDistFlory", "_bond2DArray", "_angle2DArray",
                 "_filenameAngleDist", "_dihedral2DArray", "_filenameDihedralDist",
                 "_improper2DArray", "_filenameImproperDist"]

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
        self._angle2DArray = None
        self._dihedral2DArray = None
        self._improper2DArray = None

        #Setup histograms
        self._setupHistograms()


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
    def calculate(self, begin=0, unwrap_pbc=True, type=None, ndx_filename=None, dist_label=None):

        """
         Calculate bonded distributions

        Returns:

        """

        nframes = self._trajectory.get_numframes()
        natoms = self._trajectory.natoms

        # Start calculations for each frame
        #nframes = self._trajectory.nframes
        ini = begin
        m = "\tUnwrap PBC coordinates: {}".format(unwrap_pbc)
        print(m) if self._logger is None else self._logger.info(m)

        # Main loop of frames ====================================================
        s = datetime.datetime.now()
        idx_f = 0
        X1 = np.zeros([natoms], float)
        Y1 = np.zeros([natoms], float)
        Z1 = np.zeros([natoms], float)

        for iframe in range(ini, nframes, self._stride):

            # Estimated time (Use the 10 first frames to estimate the time)
            if idx_f == 0:
                f = datetime.datetime.now()
            if idx_f == 10:
                elapsed_time = datetime.datetime.now() - f
                k = int(((nframes - ini)/self._stride))/10
                estimated_time = k*elapsed_time.total_seconds()
                m = "\tESTIMATED TIME using {0:d} frames ({1:s} seconds): {2:.2f} seconds".format \
                    (10, str(elapsed_time.total_seconds()), estimated_time)
                print(m) if self._logger is None else self._logger.info(m)

            # Write info
            if iframe%self._freq == 0:
                elapsed_time = datetime.datetime.now() - s
                m = "\tIFRAME: {0:d} in {1:s} seconds".format \
                    (iframe, str(elapsed_time.total_seconds()))
                print(m) if self._logger is None else self._logger.info(m)

            # If pbc is false it is assumed that the trajectory is unwrapped
            if unwrap_pbc:
                self._unwrap_coordinates(iframe)
            else:
                self._coords_unwrap = self._trajectory.universe.trajectory[iframe].positions

            X1[:] = self._coords_unwrap[:, 0]
            Y1[:] = self._coords_unwrap[:, 1]
            Z1[:] = self._coords_unwrap[:, 2]

            if type == "bond":
                self._calculate_bond_dist_frame(ndx_filename, dist_label, X1, Y1, Z1)
            elif type == "angle":
                self._calculate_angle_dist_frame(ndx_filename, dist_label, X1, Y1, Z1)
            elif type == "dihedral":
                self._calculate_dihedral_dist_frame(ndx_filename, dist_label, X1, Y1, Z1)
            elif type =="improper":
                self._calculate_improper_dist_frame(ndx_filename, dist_label, X1, Y1, Z1)
            iframe += self._stride
            idx_f += 1

        self.writeDist(type=type, ikey=dist_label)

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

        return dihedral2DList, typedihedral_dict

    # #######################################################################
    def _calculate_dihedral_dist_frame(self, ndx_filename, dihdist_label, X1, Y1, Z1):

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

                iserror3 = self.dihDist(self._dihedral2DArray, X1, Y1, Z1)

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
    def dihDist(self,dih2DArray,X1,Y1,Z1):

        iserror = 1
        if len(dih2DArray) != 0:
            iserror = pag.dihDistC(dih2DArray, X1, Y1, Z1, self._dihHist)
        return iserror

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
            with open(self._filenameDihedralDist, 'w') as f:
                for i in range(0,self._maxbinDih):
                    if np.sum(self._dihHist) != 0:
                        f.writelines(str(self._ddist[i])+" "+str(self._dihHist[i])+" "+\
                                     str(self._dihHist[i]*(1.0/np.sum(self._dihHist))*(1.0/self._deltaDih))+"\n")
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



        # self._fileBondDist.close()
        # self._fileAngleDist.close()
        # self._fileDihDist.close()
        # self._fileDihDistFlory.close()
        # self._fileImpDist.close()

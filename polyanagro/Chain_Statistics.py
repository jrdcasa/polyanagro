import numpy as np
import datetime
import os
import sys
import h5py
import polyanagro as pag
import topology as top
from collections import defaultdict


# ===============================================================
class Chain_Statistics(pag.Calculations):

    __slots__ = ["_removetmpfiles_ee", "_removetmpfiles_rg", "_rg2_frame", "_ree2_frame", "_ree2rg2_frame",
                 "_uree_frame", "_rEE_ACF", "_log", "_ree_max", "_ree_min", "_rg_max", "_rg_min",
                 "_molecularweigth_avg", "_all_bonds"]
    # #########################################################################
    def __init__(self, trj, dt=1, stride=1,removetmpfiles_ee=True,
                 removetmpfiles_rg=True, log=None):

        """

        ``Parameters``:
            * ``trj``:
            * ``dt``:
            * ``stride``:
            * ``removetmpfiles_ee``:
            * ``removetmpfiles_rg``:
            * ``log``:

        ``Attributes``:
            * ``self._rg2_frame`` (dict) : Radius of gyration data {Frame # :[<Rg2> stdRg2],...} in Angstroms
            * ``self._ree2_frame`` (dict) : End to end data {Frame # :[<Ree2> stdRee2],...} in Angstroms
            * ``self._ree2rg2_frame`` (dict) : <Ree2>/<Rg2> {Frame # :[<Ree2/Rg2> std<Ree2/Rg2>,...} in Angstroms
            * ``self._molecularweigth_avg`` (float) : Average Molecular weigth of the polymer chains (g/mol)


        """

        super().__init__(trj, dt=dt, stride=stride, logger=log)
        self._removetmpfiles_ee = removetmpfiles_ee
        self._removetmpfiles_rg = removetmpfiles_rg
        # Mean-square Radius of gyration for each frame (A2) --> a[iframe] = [avg, std]
        self._rg2_frame = {}
        # Mean-square end to end distance for each frame (A2) --> a[iframe] = [avg, std]
        self._ree2_frame = {}
        # Mean-square end to end distance/Mean-square Radius of guration for each frame --> a[iframe] = [avg, std]
        self._ree2rg2_frame = {}
        # Average molecular weight per chain
        self._molecularweigth_avg = 1.0
        self._uree_frame = None
        self._rEE_ACF = None
        self._all_bonds = None

        self._ree_max = sys.float_info.min
        self._ree_min = sys.float_info.max
        self._rg_max = sys.float_info.min
        self._rg_min = sys.float_info.max

    # #########################################################################
    def calculate(self, listendtoend, isree=True, isrg=True, iscn=True, begin=0,
                  molecularweight=True, unwrap_pbc=True, diroutput="./",
                  acfE2E=False, distributions=False,
                  backbone_list_atoms=None, isbbatom=None, calc_Cn_bonds_distances=False,
                  single_Cn_unitvector=False):

        """
        Calculate dimensions of the polymer chains or molecules

            * Molecular weight calculation: It is supposed that the molecular Mw is constant in the whole trajectory

        ``Parameters``:
            * ``write_rg``:
            * ``molecularweight``(boolean): Calculate the molecular average weight of the polymer
            * ``pbc``:
            * ``diroutput``:
            * ``listendtoend``:
            * ``acfE2E``:
            * ``distributions``:
            * ``begin``: First frame in the trajectory to be used

        ``Return``:
            * ``None``
        """

        nchains = len(self._nmols_array)
        nframes = self._trajectory.get_numframes()

        # End to end cannot be calculated without listendtoend
        if isree and listendtoend is None:
            m = "\tEnd to end calculation is deactivated, listendtoend is None\n"
            print(m) if self._logger is None else self._logger.info(m)
            isree = False
        # Cn cannot be calculated without branch and backbone information
        if (iscn and isbbatom is None) or\
            (iscn and backbone_list_atoms is None):
            m = "\tCn calculation is deactivated, either not info about branch or backbone\n"
            print(m) if self._logger is None else self._logger.info(m)
            iscn = False

        if diroutput[-1] != "/":
            diroutput = diroutput+"/"

        # Write messages
        start_time = datetime.datetime.now()
        if molecularweight:
            m = ""
            m += "\t*** Molecular weigth...\n"
            print(m) if self._logger is None else self._logger.info(m)
            mass_by_chain = np.zeros(nchains)
            for ich in range(nchains):
                m = 0.0
                for el in self._nmols_array[ich]:
                    m+=self._trajectory.topology.mass[el]
                mass_by_chain[ich] = m
            self._molecularweigth_avg = np.mean(mass_by_chain)
            m = "\tAverage Molecular Weight (g/mol) : {0:.2f}\n".format(self._molecularweigth_avg)
            print(m) if self._logger is None else self._logger.info(m)

        # log messages
        m = ""
        m += "\t*** Calculating Radius of gyration...\n"
        print(m) if self._logger is None else self._logger.info(m)
        m = "\tRadius of gyration to be written to {0:s} \n".format(diroutput+"Rg.dat")
        print(m) if self._logger is None else self._logger.info(m)
        # log messages
        m  = ""
        m += "\t*** Calculating End to End distance...\n\n"
        m += "\tEnd to End distance to be written to {0:s} \n".format(diroutput+"Ree.dat")
        print(m) if self._logger is None else self._logger.info(m)

        # Start calculations for each frame
        #nframes = self._trajectory.nframes
        ini = begin
        m = "\tUnwrap PBC coordinates: {}".format(unwrap_pbc)
        print(m) if self._logger is None else self._logger.info(m)

        # Check if Cn with unit vectors can be used
        # if self._trajectory.topology._topologyfile.split(".")[-1] != "psf":
        #     m = "\tCalculation of Cn with unit vectors have been required but\n"
        #     m += "\ttopology is not psf.\n"
        #     print(m) if self._logger is None else self._logger.info(m)
        #     single_Cn_unitvector = False

        # For all frames
        if iscn:
            self._all_bonds = self._trajectory.topology.get_allbonds()

        for iframe in range(ini, nframes, self._stride):

            # Write info
            if iframe%self._freq == 0:
                elapsed_time = datetime.datetime.now() - start_time
                m = "\tIFRAME: {0:d} in {1:s} seconds".format \
                    (iframe, str(elapsed_time.total_seconds()))
                print(m) if self._logger is None else self._logger.info(m)

            # If pbc is false it is assumed that the trajectory is unwrapped
            if unwrap_pbc:
                self._unwrap_coordinates(iframe)
            else:
                self._coords_unwrap = self._trajectory.universe.trajectory[iframe].positions

            if isrg:
                self._single_rg(iframe, filename=diroutput+"Rg.dat", tmpfileindividualchain=distributions)
            if isree:
                self._single_ree(iframe, filename=diroutput+"Ree.dat",
                             listendtoend=listendtoend,tmpfileindiviualchain=distributions)
            if isrg and isree:
                self._single_ree2rg2(iframe, filename=diroutput + "Ree2Rg2.dat")
            if iscn:
                self._single_Cn(iframe, backbone_list_atoms, isbbatom,
                            filename=diroutput + "Cn.dat", calc_distances=calc_Cn_bonds_distances,
                            single_Cn_unitvector=single_Cn_unitvector)

            iframe += self._stride

         # Write final messages and time
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\n\tTIME(Chain dimensions): {0:s} seconds\n".\
             format(str(elapsed_time.total_seconds()))
        m += "\t*** End Calculating Chain Dimensions.\n"
        print(m) if self._logger is None else self._logger.info(m)

        if acfE2E:
            # Start ndarray to store the end to end vectors for autocorrelation function
            self._uree_frame = np.zeros([3, nchains, nframes], dtype=np.float32)
            self._rEE_ACF = np.zeros([nframes], dtype=np.float32)

            # Write messages
            start_time = datetime.datetime.now()

            m = "\tCalculating End-to-End autocorrelation vector..."
            print(m) if self._logger is None else self._logger.info(m)

            self._acf_e2e(filename=diroutput+"rEE_ACF.dat")

            # Write final messages and time
            end_time = datetime.datetime.now()
            elapsed_time = end_time - start_time
            m = "\tTIME(End-to-End autocorrelation vector): {0:s} seconds\n".format(str(elapsed_time.total_seconds()))
            m += "\t*** End Calculating End-to-End autocorrelation vector.\n"
            print(m) if self._logger is None else self._logger.info(m)


    #     if distributions:
    #         with open(".tmp_rg_distances.dat", 'a') as f:
    #             f.writelines("Rg_min_max: {0:.4f} {1:.4f}\n".format(self._rg_min, self._rg_max))
    #         if not listendtoend is None:
    #             with open(".tmp_ree_distances.dat", 'a') as f:
    #                 f.writelines("Ree_min_max: {0:.4f} {1:.4f}\n".format(self._ree_min, self._ree_max))
    #
    #         self._distribution(".tmp_ree_distances.dat")
    #         self._distribution(".tmp_rg_distances.dat")
    #
    #         print("Calculate distributions TO BE DONE!!!!!!")
    #


    # #########################################################################
    def _single_rg(self, iframe, filename="Rg.dat", write_result=True,
                   ismassweight=False, tmpfileindividualchain=False):

        mass = np.array(self._trajectory.topology.mass, dtype=np.float32)
        nchains = len(self._nmols_array)
        rgsq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)

        if ismassweight:
            pag.calc_rg_openmp_massweigth(self._nmols_array, mass, self._coords_unwrap, rgsq_ich_iframe)
        else:
            pag.calc_rg_openmp(self._nmols_array, mass, self._coords_unwrap, rgsq_ich_iframe)

        # Calculate the average
        rg2_avg = 0.0
        rg2_std = 0.0

        # Calculate the avaerage
        rg2_avg = np.mean(rgsq_ich_iframe[:,3])
        # Calculate the standard deviation
        rg2_std = np.std(rgsq_ich_iframe[:,3])

        self._rg2_frame[iframe] = [rg2_avg, rg2_std]

        # Write a tmp file containing the radius of gyration for each chain.
        # Each line corresponds to all chains in the current frame
        # This is useful to calculate the radius of gyration distribution
        if tmpfileindividualchain:
            try:
                if self._removetmpfiles_rg:
                    os.remove(".tmp_rg_distances.dat")
            except FileNotFoundError:
                pass
            with open(".tmp_rg_distances.dat",'a') as f:
                self._removetmpfiles_rg = False
                line = "{} ".format(iframe)
                for ich in range(nchains):
                    rg = rgsq_ich_iframe[ich,3]**0.5
                    line += " {0:.4f} ".format(rg)
                    self._rg_max = rg if rg > self._rg_max else self._rg_max
                    self._rg_min = rg if rg < self._rg_min else self._rg_min
                line += "\n"
                f.writelines(line)

        # ====== DATA FOR Rg =======
        if write_result:
            i = iframe * self._dt
            if iframe == 0:
                with open(filename,'w') as f:
                    f.writelines("#  iFrame          Time(ps)     <Rg^2>(A^2)     stdRg^2(A^2)     <Rg^2>/M(molnm2/kg)\n")
                    f.writelines("#===================================================================================\n")

                    f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}   {4:>10.3f}\n".format(
                                                    iframe,i,
                                                    self._rg2_frame[iframe][0],
                                                    self._rg2_frame[iframe][1],
                                                    self._rg2_frame[iframe][0]*10/self._molecularweigth_avg))
            else:
                with open(filename,'a') as f:
                    f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}   {4:>10.3f}\n".format(
                                                    iframe, i,
                                                    self._rg2_frame[iframe][0],
                                                    self._rg2_frame[iframe][1],
                                                    self._rg2_frame[iframe][0]*10/self._molecularweigth_avg))

        return rg2_avg, rg2_std

    # #########################################################################
    def _single_ree(self, iframe, filename="Ree.dat", write_result=True,
                    listendtoend=None, tmpfileindiviualchain=False):

        nchains = len(self._nmols_array)
        reesq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)
        ref_head = np.zeros((nchains,3))
        conf_tail = np.zeros((nchains,3))

        # Calculate the average
        ree2_avg = 0.0
        ree2_std = 0.0

        for item in listendtoend:
            ich, ihead, itail = item
            ref_head[ich,:] = self._coords_unwrap[ihead,:]
            conf_tail[ich,:] = self._coords_unwrap[itail,:]

        dij, reex, reey, reez = top.distance_diagonal_array(ref_head, conf_tail, openmp=True)

        # Calculate the average
        reesq_ich_iframe[:,3] = dij[:]*dij[:]
        ree2_avg = np.mean(reesq_ich_iframe[:,3])
        # Calculate the standard deviation
        ree2_std = np.std(reesq_ich_iframe[:,3])

        for ich in range(nchains):
            if not self._uree_frame is None:
                self._uree_frame[0, ich, iframe] = reex[ich] / dij[ich]
                self._uree_frame[1, ich, iframe] = reey[ich] / dij[ich]
                self._uree_frame[2, ich, iframe] = reez[ich] / dij[ich]

        self._ree2_frame[iframe] = [ree2_avg, ree2_std]

        # ====== DATA FOR Ree DISTRIBUTIONS =======
        # Write a tmp file containing all end-to-end distances for each chain.
        # Each line corresponds to all chains in the current frame
        # This is useful to calculate the end-to-end distribution
        if tmpfileindiviualchain:
            try:
                if self._removetmpfiles_ee:
                    os.remove(".tmp_ree_distances.hdf5")
            except FileNotFoundError:
                pass
            # with open(".tmp_ree_distances.dat",'a') as f:
            #     self._removetmpfiles_ee = False
            #     line = "{} ".format(iframe_total)
            #     for ich_ee in dij:
            #         line += " {0:.4f} ".format(ich_ee)
            #         self._ree_max = ich_ee if ich_ee > self._ree_max else self._ree_max
            #         self._ree_min = ich_ee if ich_ee < self._ree_min else self._ree_min
            #     line += "\n"
            #     f.writelines(line)
            with h5py.File(".tmp_ree_distances.hdf5",'a') as f:
                self._removetmpfiles_ee = False
                arr = np.zeros(nchains)
                for ich in range(nchains):
                    arr[ich] = dij[ich]
                dset = f.create_dataset('Frame_{0:07d}'.format(iframe), data=arr)


        # ====== DATA FOR Ree =======
        if write_result:
            i = iframe * self._dt
            if iframe == 0:
                with open(filename,'w') as f:
                    f.writelines("#  iFrame          Time(ps)     <Ree^2>(A^2)     stdRee^2(A^2)     <Ree^2>/M(molnm2/kg)\n")
                    f.writelines("#======================================================================================\n")

                    f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}   {4:>10.3f}\n".format(
                                                    iframe,i,
                                                    self._ree2_frame[iframe][0],
                                                    self._ree2_frame[iframe][1],
                                                    self._ree2_frame[iframe][0]*10/self._molecularweigth_avg))
            else:
                with open(filename,'a') as f:
                    f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}   {4:>10.3f}\n".format(
                                                    iframe, i,
                                                    self._ree2_frame[iframe][0],
                                                    self._ree2_frame[iframe][1],
                                                    self._ree2_frame[iframe][0]*10/ self._molecularweigth_avg))

        return ree2_avg, ree2_std

    # #########################################################################
    def _single_ree2rg2(self, iframe, filename="Rg2Ree2.dat"):

        mass = np.array(self._trajectory.topology.mass, dtype=np.float32)
        nchains = len(self._nmols_array)
        rgsq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)

        a1 = self._ree2_frame[iframe][0] / self._rg2_frame[iframe][0]
        a2 = a1 * np.sqrt((self._ree2_frame[iframe][1] / self._ree2_frame[iframe][0]) ** 2 +
                          (self._rg2_frame[iframe][1] / self._rg2_frame[iframe][0]) ** 2)
        self._ree2rg2_frame[iframe] = (a1, a2)

        i = iframe * self._dt
        if iframe == 0:
            with open(filename,'w') as f:
                f.writelines("#  iFrame          Time(ps)     <Ree^2>/<Rg^2>  std<Ree^2>/<Rg^2>\n")
                f.writelines("#================================================================\n")

                f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f} \n".format(
                    iframe,i,
                    self._ree2rg2_frame[iframe][0],
                    self._ree2rg2_frame[iframe][1]))
        else:
            with open(filename,'a') as f:
                f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f} \n".format(
                    iframe, i,
                    self._ree2rg2_frame[iframe][0],
                    self._ree2rg2_frame[iframe][1]))

        return a1, a2

    # #########################################################################
    def _single_Cn(self, iframe,  backbone_atoms, is_bb_atoms, lavg=1.54,
                   calc_distances=False, filename="Cn.dat", single_Cn_unitvector=False):


        nchains = len(self._nmols_array)

        if backbone_atoms is None or is_bb_atoms is None:
            return None
        else:
            nat_bb = []
            natoms = len(is_bb_atoms)
            for ich in range(nchains):
                nat_bb.append(len(backbone_atoms[ich]))
            natbb_avg = np.mean(nat_bb)
            natbb_bonds_avg = natbb_avg - 1

        if not calc_distances:
                lavg2 = lavg*lavg
        else:
            #Calc all bond distances for each chain and average
            dsum = 0.0
            nelem = 0
            for ich in backbone_atoms:
                natbbch = len(ich)
                ref = np.zeros((natbbch-1,3), dtype=np.float64)
                conf = np.zeros((natbbch-1,3), dtype=np.float64)
                i = 0
                for iat in ich:
                    for ineigh in self._l_neigh_array[iat]:
                        if not is_bb_atoms[ineigh]: continue
                        if iat > ineigh: continue
                        ref[i] = self._coords_unwrap[iat,:]
                        conf[i] = self._coords_unwrap[ineigh,:]
                        i+= 1
                d = top.distance_diagonal_array(ref, conf, openmp=True)[0]
                dsum += np.sum(d)
                nelem += len(d)
            lavg = dsum/nelem
            lavg2 = lavg*lavg

        if single_Cn_unitvector:
            nbonds_bb_perch = defaultdict(int)
            all_bb_bonds = []
            for item in self._all_bonds:
                at1 = item[0]
                at2 = item[1]
                ich1 = self._trajectory.topology._iatch[at1]
                ich2 = self._trajectory.topology._iatch[at2]
                if not  is_bb_atoms[at1]: continue
                if not  is_bb_atoms[at2]: continue
                nbonds_bb_perch[ich1] += 1
                all_bb_bonds.append([at1, at2])
            all_bb_bonds = np.array(all_bb_bonds, dtype=np.int32)
            iatch = self._trajectory.topology._iatch
            nbonds_bb_max = max(nbonds_bb_perch.values())
            uux = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
            uuy = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
            uuz = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
            cn_uu = pag.unit_bond_vectors(nchains, all_bb_bonds, self._coords_unwrap, iatch, uux, uuy, uuz)
        else:
            cn_uu = 0.0

        # Write results
        i = iframe * self._dt
        if iframe == 0:
            with open(filename, 'w') as f:
                f.writelines("#  Number of bacbkone atoms:{}, calc_distances:{}, lavg :{}\n".
                             format(natbb_avg, calc_distances, lavg))
                f.writelines("#  iFrame          Time(ps)     <Ree^2>/nl^2     6.0*<Rg^2>/nl^2      Cn_formula\n")
                f.writelines("#===============================================================================\n")

                f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}    {4:>10.2f}\n".format(
                    iframe, i,
                    self._ree2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    6.0 * self._rg2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    cn_uu))
        else:
            with open(filename, 'a') as f:
                f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}    {4:>10.2f}\n".format(
                    iframe, i,
                    self._ree2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    6.0 * self._rg2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    cn_uu))

    # #########################################################################
    def odf_intra_cython(self, initframe=0, endframe=None, stride=1, filename="odf_intra.dat"):

        """
        Calculate intra chain orientation correlation (for persistence length)

        Calculate orientational correlation functions
        for intra-chain vectors as function of chemical distance.
        It also calculates components and 4th moments.
        The input is supposed to contain a multiple of nvec vectors.

        """

        nmols_array, neigbors = self._trajectory.topology.get_array_mols_neigh()
        nchains = len(nmols_array)
        all_bb_bonds, nbonds_bb_perch = self._trajectory.topology.get_all_bb_bonds()
        try:
            nbonds_bb_max = max(nbonds_bb_perch.values())
        except:
            nbonds_bb_max = 1
        coords_t0_wrapped = self._trajectory.universe.trajectory[0].positions
        box_dimensions = self._trajectory.universe.trajectory[0].dimensions[0:3]
        iatch = self._trajectory.topology._iatch
        nframes = self._trajectory.get_numframes()
        if endframe is None:
            endframe = nframes

        setup_odf_intra(nbonds_bb_max)
        uux = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
        uuy = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
        uuz = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
        for iframe in range(initframe ,endframe, stride):
            if iframe %100 == 0:
                print("{} of {}".format(iframe, nframes ))
            coords_t0_wrapped = self._trajectory.universe.trajectory[iframe].positions
            box_dimensions = self._trajectory.universe.trajectory[iframe].dimensions[0:3]
            coords_unwrap = pag.unwrap(coords_t0_wrapped, nmols_array, neigbors,
                                       box_dimensions, iframe=iframe)

            odf_intra(iframe, nchains, all_bb_bonds, coords_unwrap, iatch, uux, uuy, uuz)

        avg_write_odf_intra(nframes, nbonds_bb_max, filename)

    # # #########################################################################
    def _acf_e2e(self, filename="rEE_ACF.dat", write_result=True):

        ndumps = self._uree_frame.shape[2]
        pag.calc_acf_ete(self._uree_frame, self._rEE_ACF)

        self._dt = 1.0
        if write_result:
            for iframe_total in range(0, self._rEE_ACF.shape[0]):
                i = iframe_total*self._dt
                if iframe_total == 0:
                    with open(filename,'w') as f:
                        f.writelines("#   Time(ps)     rEE_ACF\n")
                        f.writelines("#===================================\n")
                        f.writelines("{0:>10.2f}     {1:<10.8f}   \n".format(i,
                                                        self._rEE_ACF[iframe_total]))
                else:
                    with open(filename,'a') as f:
                        f.writelines("{0:>10.2f}     {1:<10.8f} \n".format(i,
                                                        self._rEE_ACF[iframe_total]))

    # # #########################################################################
    # def _distribution(self, filename):
    #
    #     print(filename)
    #     with open(filename, 'r') as f:
    #         print(f)
    #



import numpy as np
import datetime
import os
import math
import sys
import h5py
import polyanagro as pag
import topology as top
from collections import defaultdict
import matplotlib.pyplot as plt


# ===============================================================
class Chain_Statistics(pag.Calculations):

    __slots__ = ["_removetmpfiles_ee", "_removetmpfiles_rg", "_rg2_frame", "_ree2_frame", "_ree2rg2_frame",
                 "_uree_frame", "_ree_frame", "_rEE_ACF", "_rEE2_ACF", "_log", "_ree_max", "_ree_min",
                 "_rg_max", "_rg_min",
                 "_molecularweigth_avg", "_all_bonds", "_fdist_h5py", "_nbonds_bb_max", "_all_bb_bonds",
                 "_iatch", "_cbb_avg", "_rsq_intdist", "_rsqcount_intdist", "_rsqavg_intdist", "_lavg2"]

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
        self._uree_frame = None     # Unit vector
        self._ree_frame = None      # Ree not unit vector
        self._rEE_ACF = None
        self._rEE2_ACF = None
        self._lavg2 = 0.0

        self._all_bonds = None
        self._nbonds_bb_max = None
        self._all_bb_bonds = None
        self._iatch = None
        self._cbb_avg = None

        self._ree_max = sys.float_info.min
        self._ree_min = sys.float_info.max
        self._rg_max = sys.float_info.min
        self._rg_min = sys.float_info.max

        self._fdist_h5py = None  #File in format h5py to store data to calculate distributions.


    # #########################################################################
    def calculate(self, listendtoend, isree=True, isrg=True, iscn=True, begin=0,
                  molecularweight=True, unwrap_pbc=True, diroutput="./",
                  acfE2E=False, distributions=False, isbondorientation=False,
                  backbone_list_atoms=None, isbbatom=None, calc_Cn_bonds_distances=False,
                  single_Cn_unitvector=False, isinternalchaindist=True, isrgmass=False, isodf=False):

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
        if self._stride == 1:
            nframes_analysed = nframes
        else:
            nframes_analysed = int(((nframes - begin) / self._stride)) + 1
        m = "\t Num of frames to analyse: {}".format(nframes_analysed)
        print(m) if self._logger is None else self._logger.info(m)

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

        # ODF cannot be calculated without branch and backbone information
        if (isodf and isbbatom is None) or\
           (isodf and backbone_list_atoms is None):
            m = "\tODF calculation is deactivated, either not info about branch or backbone\n"
            print(m) if self._logger is None else self._logger.info(m)
            isodf = False

        # Internal chain distances cannot be calculated without branch and backbone information
        if isinternalchaindist and backbone_list_atoms is None:
            m = "\tInternal chain distances are deactivated, not info about backbone\n"
            print(m) if self._logger is None else self._logger.info(m)
            isinternalchaindist = False

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
            m = "\tAverage Molecular Weigth (g/mol) : {0:.2f}\n".format(self._molecularweigth_avg)
            print(m) if self._logger is None else self._logger.info(m)

        # log messages
        m = ""
        m += "\t*** Calculating Radius of gyration..."
        print(m) if self._logger is None else self._logger.info(m)
        m = "\tRadius of gyration to be written to {0:s} \n".format(diroutput+"Rg.dat")
        print(m) if self._logger is None else self._logger.info(m)
        # log messages
        m  = ""
        m += "\t*** Calculating End to End distance...\n"
        m += "\tEnd to End distance to be written to {0:s} \n".format(diroutput+"Ree.dat")
        print(m) if self._logger is None else self._logger.info(m)

        # Start calculations for each frame
        #nframes = self._trajectory.nframes
        ini = begin
        m = "\tUnwrap PBC coordinates: {}".format(unwrap_pbc)
        print(m) if self._logger is None else self._logger.info(m)

        # For all frames
        if (iscn or isbondorientation or isinternalchaindist or isodf) and isree and len(listendtoend) != 0:

            s = datetime.datetime.now()
            nbonds_bb_perch = defaultdict(int)
            all_bb_bonds = []
            # For each chain walk from the head to end of
            # the chain through the backbone atoms.
            # Bonds in the backbon have to be ordered from head to tail
            isvisited = [False]*self._trajectory.topology.natoms
            ichatbb = np.zeros(self._trajectory.topology.natoms, dtype=np.int32)
            ichatbb[:] = -1
            for item in listendtoend:
                ich = item[0]
                ihead = item[1]
                queue = [ihead]
                while queue:
                    iat = queue.pop()
                    isvisited[iat] = True
                    neigh = self._trajectory.topology.get_neighbours(iat)
                    for jat in neigh:
                        if isbbatom[jat] and not isvisited[jat]:
                            ich1 = self._trajectory.topology._iatch[iat]
                            ich2 = self._trajectory.topology._iatch[jat]
                            nbonds_bb_perch[ich1] += 1
                            all_bb_bonds.append([iat, jat])
                            ichatbb[iat] = jat
                            iat = jat
                            queue.append(iat)
                            break
            self._all_bb_bonds = np.array(all_bb_bonds, dtype=np.int32)
            self._iatch = self._trajectory.topology._iatch
            self._nbonds_bb_max = max(nbonds_bb_perch.values())
            #cJ self._cbb_avg = np.zeros([nframes, self._nbonds_bb_max])
            self._cbb_avg = np.zeros([nframes_analysed, self._nbonds_bb_max])

            elapsed_time = datetime.datetime.now() - start_time
            m = "\tPreparing backbone information in {0:10.3f} seconds".format(elapsed_time.total_seconds())
            print(m) if self._logger is None else self._logger.info(m)

        else:
            m = "\tCn, bond orientation, odf and internal distances calculations " \
                "are deactivated, either not info about branch or backbone\n"
            print(m) if self._logger is None else self._logger.info(m)
            iscn = False
            isbondorientation = False
            isinternalchaindist = False
            isodf = False

        # Delete temporal files for distributions if they exist and create new datasets
        if distributions:
            s = datetime.datetime.now()
            try:
                   os.remove(".tmp_distributions.hdf5")
            except FileNotFoundError:
                pass
            self._fdist_h5py = h5py.File('.tmp_distributions.hdf5', 'w')
            #cJ self._fdist_h5py.create_dataset("rg", (nframes, nchains))
            #cJ self._fdist_h5py.create_dataset("ree", (nframes, nchains))
            self._fdist_h5py.create_dataset("rg", (nframes_analysed, nchains))
            self._fdist_h5py.create_dataset("ree", (nframes_analysed, nchains))


            elapsed_time = datetime.datetime.now() - s
            m = "\tPreparing distribution information in {0:10.3f} seconds".format(elapsed_time.total_seconds())
            print(m) if self._logger is None else self._logger.info(m)


        if acfE2E:
            #cJ self._uree_frame = np.zeros([3, nchains, nframes], dtype=np.float32)
            self._uree_frame = np.zeros([3, nchains, nframes_analysed], dtype=np.float32)
            self._ree_frame = np.zeros([3, nchains, nframes_analysed], dtype=np.float32)

        if isodf:
            pag.setup_odf_intra(self._nbonds_bb_max)
            uux = np.zeros((nchains, self._nbonds_bb_max), dtype=np.float32)
            uuy = np.zeros((nchains, self._nbonds_bb_max), dtype=np.float32)
            uuz = np.zeros((nchains, self._nbonds_bb_max), dtype=np.float32)

        if isinternalchaindist:

            s = datetime.datetime.now()
            # Prepare internal chain distances
            nchains = len(self._nmols_array)
            head_array = np.zeros(nchains, dtype=np.int32)
            for ich in range(0,nchains):
                head_array[ich] = listendtoend[ich][1]
            self._rsq_intdist = np.zeros(self._nbonds_bb_max + 1, dtype=np.float64)
            self._rsqcount_intdist = np.zeros(self._nbonds_bb_max + 1, dtype=np.int32)
            self._rsqavg_intdist = np.zeros(self._nbonds_bb_max + 1, dtype=np.float64)
            pag.setup_internal_distances(self._nbonds_bb_max)

            elapsed_time = datetime.datetime.now() - s
            m = "\tPreparing internal chain distances information in {0:10.3f} seconds".format(elapsed_time.total_seconds())
            print(m) if self._logger is None else self._logger.info(m)


        # Main loop of frames ====================================================
        s = datetime.datetime.now()
        idx_f = 0
        for iframe in range(ini, nframes, self._stride):

            # Estimated time (Use the 100 first frames to estimate the time)
            if idx_f == 0:
                f = datetime.datetime.now()
            if idx_f == 100:
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
            if unwrap_pbc:
                self._unwrap_coordinates(iframe)
            else:
                self._coords_unwrap = self._trajectory.universe.trajectory[iframe].positions

            if isrg:
                fname = os.path.join(diroutput, "Rg.dat")
                self._single_rg(iframe, idx_f, filename=fname, tmpfileindividualchain=distributions)

            if isrgmass:
                self._single_rg(iframe, idx_f, filename="Rg_mass.dat", ismassweight=True)

            if isree:
                fname = os.path.join(diroutput, "Ree.dat")
                self._single_ree(iframe, idx_f, filename=fname,
                             listendtoend=listendtoend,tmpfileindividualchain=distributions)
            if isrg and isree:
                fname = os.path.join(diroutput, "Ree2Rg2.dat")
                self._single_ree2rg2(iframe, filename=fname)

            if iscn:
                fname = os.path.join(diroutput, "Cn.dat")
                self._single_Cn(iframe, backbone_list_atoms, isbbatom,
                            filename=fname, calc_distances=calc_Cn_bonds_distances,
                            single_Cn_unitvector=single_Cn_unitvector)
            if isbondorientation:
                self._single_bond_orientation(iframe, backbone_list_atoms, isbbatom)

            if isinternalchaindist:
                #self._insternalchaindist(iframe, head_array, ichatbb)
                #self._insternalchaindist_python(iframe, head_array, ichatbb)
                pag.insternalchaindist_cython(nchains, self._coords_unwrap, head_array, ichatbb, self._rsq_intdist,
                                              self._rsqcount_intdist)
            if isodf:
                pag.odf_intra(iframe, nchains, self._all_bb_bonds, self._coords_unwrap, self._iatch, uux, uuy, uuz)

            iframe += self._stride
            idx_f += 1

        if isbondorientation:
            self._write_bond_orientation()

        if isinternalchaindist:
            self._write_isinternalchaindist(os.path.join(diroutput,"cn_internal_distances.dat"))

        # Write final messages and time
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\n\tTIME(Chain dimensions): {0:s} seconds\n".\
             format(str(elapsed_time.total_seconds()))
        m += "\t*** End Calculating Chain Dimensions.\n"
        print(m) if self._logger is None else self._logger.info(m)

        if acfE2E:
            # Start ndarray to store the end to end vectors for autocorrelation function

            #cJ self._rEE_ACF = np.zeros([nframes], dtype=np.float32)
            self._rEE_ACF = np.zeros([nframes_analysed], dtype=np.float32)
            self._rEE2_ACF = np.zeros([nframes_analysed], dtype=np.float32)

            # Write messages
            start_time = datetime.datetime.now()
            m = "\t*** Calculating End-to-End autocorrelation vector... ({})".format(start_time)
            print(m) if self._logger is None else self._logger.info(m)

            self._acf_e2e(filename=diroutput+"rEE_ACF.dat")

            # Write final messages and time
            end_time = datetime.datetime.now()
            elapsed_time = end_time - start_time
            m = "\tTIME(End-to-End autocorrelation vector): {0:s} seconds\n".format(str(elapsed_time.total_seconds()))
            m += "\t*** End Calculating End-to-End autocorrelation vector.\n"
            print(m) if self._logger is None else self._logger.info(m)

        if distributions:
            start_time = datetime.datetime.now()
            m = "\t*** Calculating distributions..."
            print(m) if self._logger is None else self._logger.info(m)

            self._distribution(".tmp_distributions.hdf5", binsize=1, plothistogram=False)
            os.remove(".tmp_distributions.hdf5")

            end_time = datetime.datetime.now()
            elapsed_time = end_time - start_time
            m = "\tTIME(Distributions): {0:s} seconds\n".format(str(elapsed_time.total_seconds()))
            m += "\t*** End distributions calculation.\n"
            print(m) if self._logger is None else self._logger.info(m)

        if isodf:
            fname = os.path.join(diroutput, "odf_intra.dat")
            #pag.avg_write_odf_intra(nframes, self._nbonds_bb_max, fname)
            pag.avg_write_odf_intra(nframes_analysed, self._nbonds_bb_max, fname)

    # #########################################################################
    def statistics(self, fraction_trj, diroutput="./"):

        # Rg avg =========================================
        fname = os.path.join(diroutput, "Rg.dat")
        data_array = self._extract_data_for_avg(fname, cols=[0, 2, 3])
        Rg_avg, Rg_std, Rg_avgstd = self._calc_avg(data_array, fraction_trj)

        # Ree avg =========================================
        fname = os.path.join(diroutput, "Ree.dat")
        if os.path.isfile(fname):
            data_array = self._extract_data_for_avg(fname, cols=[0, 2, 3])
            Ree_avg, Ree_std, Ree_avgstd = self._calc_avg(data_array, fraction_trj)
            # Ree2Rg avg =========================================
            fname = os.path.join(diroutput, "Ree2Rg2.dat")
            data_array = self._extract_data_for_avg(fname, cols=[0, 2, 3])
            ReeRg_avg, ReeRg_std, ReeRg_avgstd = self._calc_avg(data_array, fraction_trj)
        else:
            Ree_avg = -999999.99
            Ree_avgstd = -999999.99
            ReeRg_avg = -999999.99
            ReeRg_avgstd = -999999.99

        # Cn avg =========================================
        fname = os.path.join(diroutput, "Cn.dat")
        if os.path.isfile(fname):
            data_array = self._extract_data_for_avg(fname, cols=[0, 2, None])
            Cn_avg, Cn_std, Cn_avgstd = self._calc_avg(data_array, fraction_trj)
        else:
            Cn_avg = -999999.99
            Cn_avgstd = -999999.99

        m = "\t*** Average values ***\n"
        m += "\t\t Rg^2_avg     = {0:10.2f} +- {1:10.2f} (angstroms)\n".format(Rg_avg, Rg_avgstd)
        m += "\t\t Ree^2_avg    = {0:10.2f} +- {1:10.2f} (angstroms)\n".format(Ree_avg, Ree_avgstd)
        m += "\t\t Ree^2/Rg^2   = {0:10.2f} +- {1:10.2f} (angstroms)\n".format(ReeRg_avg, ReeRg_avgstd)
        m += "\t\t Cn_avg       = {0:10.2f} +- {1:10.2f} \n".format(Cn_avg, Cn_avgstd)
        print(m) if self._logger is None else self._logger.info(m)

        # Template plots dimensions
        dict_avg = defaultdict()
        dict_avg["Rg"] = [Rg_avg, Rg_avgstd]
        dict_avg["Ree"] = [Ree_avg, Ree_avgstd]
        dict_avg["Ree2Rg2"] = [ReeRg_avg, ReeRg_avgstd]
        dict_avg["Cn"] = [Cn_avg, Cn_avgstd]
        fnamegnu = os.path.join(diroutput, "gnuplot_dimensions.gnu")
        self._gnuplot_template_dimensions(fnamegnu, dict_avg)

        # Template plots distributions
        fnamegnu = os.path.join(diroutput, "gnuplot_distributions.gnu")
        self._gnuplot_template_distributions(fnamegnu, dict_avg)

        # Template plots characteristic ratio
        fnamegnu = os.path.join(diroutput, "gnuplot_charratio.gnu")
        self._gnuplot_template_charratio(fnamegnu, dict_avg)

    # #########################################################################
    def _single_rg(self, iframe, idx_fr, filename="Rg.dat", write_result=True,
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
            for ich in range(nchains):
                rg = rgsq_ich_iframe[ich, 3] ** 0.5
                self._fdist_h5py['rg'][idx_fr, ich] = rg
                self._rg_max = rg if rg > self._rg_max else self._rg_max
                self._rg_min = rg if rg < self._rg_min else self._rg_min


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
    def _single_ree(self, iframe, idx_fr, filename="Ree.dat", write_result=True,
                    listendtoend=None, tmpfileindividualchain=False):

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

        # dij [nchains] --> End-to-end distance in this frame
        # reex[nchains] --> X-coordinate of the end-to-end distance vector
        # reey[nchains] --> Y-coordinate of the end-to-end distance vector
        # reez[nchains] --> Z-coordinate of the end-to-end distance vector
        dij, reex, reey, reez = top.distance_diagonal_array(ref_head, conf_tail, openmp=True)

        # Calculate the average
        reesq_ich_iframe[:,3] = dij[:]*dij[:]
        ree2_avg = np.mean(reesq_ich_iframe[:,3])
        # Calculate the standard deviation
        ree2_std = np.std(reesq_ich_iframe[:,3])

        for ich in range(nchains):
            if not self._uree_frame is None:
                #cJ self._uree_frame[0, ich, iframe] = reex[ich] / dij[ich]
                #cJ self._uree_frame[1, ich, iframe] = reey[ich] / dij[ich]
                #cJ self._uree_frame[2, ich, iframe] = reez[ich] / dij[ich]
                self._uree_frame[0, ich, idx_fr] = reex[ich] / dij[ich]
                self._uree_frame[1, ich, idx_fr] = reey[ich] / dij[ich]
                self._uree_frame[2, ich, idx_fr] = reez[ich] / dij[ich]
                # cJ 28-01-2025
                self._ree_frame[0, ich, idx_fr] = reex[ich]
                self._ree_frame[1, ich, idx_fr] = reey[ich]
                self._ree_frame[2, ich, idx_fr] = reez[ich]

        self._ree2_frame[iframe] = [ree2_avg, ree2_std]

        # ====== DATA FOR Ree DISTRIBUTIONS =======
        # Write a tmp file containing all end-to-end distances for each chain.
        # Each line corresponds to all chains in the current frame
        # This is useful to calculate the end-to-end distribution
        if tmpfileindividualchain:
            for ich in range(nchains):
                #cJ self._fdist_h5py['ree'][iframe,ich] = dij[ich]
                self._fdist_h5py['ree'][idx_fr ,ich] = dij[ich]
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
    def _single_ree2rg2(self, iframe, filename="Rg2Ree2.dat", tmpfileindividualchain=False):

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

        self._lavg2 = lavg2

        if single_Cn_unitvector:
            uux = np.zeros((nchains, self._nbonds_bb_max), dtype=np.float32)
            uuy = np.zeros((nchains, self._nbonds_bb_max), dtype=np.float32)
            uuz = np.zeros((nchains, self._nbonds_bb_max), dtype=np.float32)
            cn_uu = pag.unit_bond_vectors(nchains, self._all_bb_bonds, self._coords_unwrap,
                                          self._iatch, uux, uuy, uuz)
        else:
            cn_uu = 0.0

        # Write results
        i = iframe * self._dt
        if iframe == 0:
            with open(filename, 'w') as f:
                f.writelines("#  Number of backbone atoms:{}, calc_distances:{}, lavg :{}\n".
                             format(natbb_avg, calc_distances, lavg))
                f.writelines("#  iFrame          Time(ps)     <Ree^2>/nl^2     6.0*<Rg^2>/nl^2      Cn_formula    n    l^2(A)\n")
                f.writelines("#==============================================================================================\n")

                f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}    {4:>10.2f}   {5:>10.1f}    {6:>6.3f}\n".format(
                    iframe, i,
                    self._ree2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    6.0 * self._rg2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    cn_uu, natbb_bonds_avg, lavg2))
        else:
            with open(filename, 'a') as f:
                f.writelines("{0:10d}   {1:>12.3f}   {2:>10.2f}   {3:>10.2f}    {4:>10.2f}   {5:>10.1f}    {6:>6.3f}\n".format(
                    iframe, i,
                    self._ree2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    6.0 * self._rg2_frame[iframe][0] / (natbb_bonds_avg * lavg2),
                    cn_uu, natbb_bonds_avg, lavg2))


    # #########################################################################
    def _acf_e2e(self, filename="rEE_ACF.dat", write_result=True):

        ndumps = self._uree_frame.shape[2]
        # Using the Doros' definiton
        pag.calc_acf_ete(self._uree_frame, self._rEE_ACF)
        # Using the defniiton in Polymer, 268, 2023, 125677 (J. Ramirez)
        pag.calc_acf2_ete(self._ree_frame, self._rEE2_ACF)

        if write_result:
            for iframe_total in range(0, self._rEE_ACF.shape[0]):
                i = iframe_total*self._dt*self._stride
                if iframe_total == 0:
                    with open(filename,'w') as f:
                        f.writelines("#   iFrame     Time(ps)     rEE_ACF\n")
                        f.writelines("#===================================\n")
                        f.writelines("{0:10d}     {1:>10.2f}     {2:<10.8f}   \n".format(iframe_total, i,
                                                        self._rEE_ACF[iframe_total]))
                else:
                    with open(filename,'a') as f:
                        f.writelines("{0:10d}     {1:>10.2f}     {2:<10.8f}   \n".format(iframe_total, i,
                                                        self._rEE_ACF[iframe_total]))

            for iframe_total in range(0, self._rEE2_ACF.shape[0]):
                i = iframe_total*self._dt*self._stride
                if iframe_total == 0:
                    with open("rEE_ACF2.dat",'w') as f:
                        f.writelines("#   iFrame     Time(ps)     rEE_ACF\n")
                        f.writelines("#===================================\n")
                        f.writelines("{0:10d}     {1:>10.2f}     {2:<10.8f}   \n".format(iframe_total, i,
                                                        self._rEE2_ACF[iframe_total]))
                else:
                    with open("rEE_ACF2.dat",'a') as f:
                        f.writelines("{0:10d}     {1:>10.2f}     {2:<10.8f}   \n".format(iframe_total, i,
                                                        self._rEE2_ACF[iframe_total]))

    # #########################################################################
    @staticmethod
    def _distribution(filename, binsize=1, plothistogram=False):

        # End to end distributions
        with h5py.File(".tmp_distributions.hdf5", 'r') as f:
            data_ree = np.array(f['ree']).flatten()
            max_value = math.ceil(data_ree.max())
            min_value = math.floor(data_ree.min())
            range_data = max_value-min_value
            nbins = int(range_data/binsize*0.5)
            hist_ree, bins = np.histogram(data_ree, bins=nbins, range=(min_value, max_value), density=True)
            newbins_ree = np.zeros(nbins)
            for ibin in range(1, len(bins)):
                newbins_ree[ibin-1] = (bins[ibin]-bins[ibin-1])/2.0+bins[ibin-1]
            # Write data to file
            with open("Ree_distribution.dat", "w") as f:
                f.writelines("# ibin(Angstroms) PDF(A-1)\n")
                f.writelines("#=========================\n")
                for ibin in range(0, len(newbins_ree)):
                    f.writelines("{0:10.4f} {1:10.4f}\n".format(newbins_ree[ibin], hist_ree[ibin]))

        # Radius of gyration distributions
        with h5py.File(".tmp_distributions.hdf5", 'r') as f:
            data_rg = np.array(f['rg']).flatten()
            max_value = math.ceil(data_rg.max())
            min_value = math.floor(data_rg.min())
            range_data = max_value-min_value
            nbins = int(range_data/binsize)
            hist_rg, bins = np.histogram(data_rg, bins=nbins, range=(min_value, max_value), density=True)
            newbins_rg = np.zeros(nbins)
            for ibin in range(1, len(bins)):
                newbins_rg[ibin-1] = (bins[ibin]-bins[ibin-1])/2.0+bins[ibin-1]
            # Write data to file
            with open("Rg_distribution.dat", "w") as f:
                f.writelines("# ibin(Angstroms) PDF(A-1)\n")
                f.writelines("#=========================\n")
                for ibin in range(0, len(newbins_rg)):
                    f.writelines("{0:10.4f} {1:10.4f}\n".format(newbins_rg[ibin], hist_rg[ibin]))

        # Radius of end-to-end^2 /Radius of gyration^2 distributions
        data_ree2rg2 = (data_ree*data_ree) / (data_rg*data_rg)
        max_value = math.ceil(data_ree2rg2.max())
        min_value = math.floor(data_ree2rg2.min())
        range_data = max_value-min_value
        nbins = int(range_data/binsize)
        hist_ree2rg2, bins = np.histogram(data_ree2rg2, bins=nbins, range=(min_value, max_value), density=True)
        newbins_ree2rg2 = np.zeros(nbins)
        for ibin in range(1, len(bins)):
            newbins_ree2rg2[ibin-1] = (bins[ibin]-bins[ibin-1])/2.0+bins[ibin-1]
        # Write data to file
        with open("Ree2Rg2_distribution.dat", "w") as f:
            f.writelines("# ibin(Angstroms) PDF(A-1)\n")
            f.writelines("#=========================\n")
            for ibin in range(0, len(newbins_ree2rg2)):
                f.writelines("{0:10.4f} {1:10.4f}\n".format(newbins_ree2rg2[ibin], hist_ree2rg2[ibin]))

        if plothistogram:
            fig, (ax1, ax2) = plt.subplots(1, 2)
            ax1.set_title('Ree distribution')
            ax1.set_xlabel(r'$R_{EE}$ ($\mathring{A}$)')
            ax1.set_ylabel(r'PDF')
            ax1.plot(newbins_ree, hist_ree)

            ax2.set_title('Rg distribution')
            ax2.set_xlabel(r'$R_g$ ($\mathring{A}$)')
            ax2.set_ylabel(r'PDF')
            ax2.plot(newbins_rg, hist_rg)
            plt.show()

    # #########################################################################
    def _single_bond_orientation(self, iframe,  backbone_atoms, is_bb_atoms):

        nchains = len(self._nmols_array)
        cbb = np.zeros(self._nbonds_bb_max, dtype=np.float32)
        pag.bond_bond_correlation(nchains, self._nbonds_bb_max,
                                  self._all_bb_bonds, self._coords_unwrap, self._iatch, cbb)

        self._cbb_avg[iframe, :] = cbb
        pass

    # #########################################################################
    def _write_bond_orientation(self):

        m = self._cbb_avg.shape[1]
        filename = "bond-bond-correlation.dat"

        with open(filename, 'w') as f:
            f.writelines("#  Number of bonds cbb  std_cbb\n")
            f.writelines("#  ============================\n")

        # Average over all frames each m value (m = j - i)
        for idx in range(1, m):
            avg = np.mean(self._cbb_avg[:,idx])
            std = np.std(self._cbb_avg[:,idx])
            with open(filename, 'a') as f:
                f.writelines("{0:8d} {1:10.4f} {2:10.4f}\n".format(idx, avg, std))

    # #########################################################################
    def _insternalchaindist(self, iframe, head_array, ichatbb):

        nchains = len(self._nmols_array)
        pag.internal_distances_iframe(nchains, head_array, ichatbb, self._coords_unwrap)

    # #########################################################################
    def _insternalchaindist_python(self, iframe, head_array, ichatbb):

        r = np.zeros(3, dtype=np.float32)

        nchains = len(self._nmols_array)
        for ich in range(0, nchains):
            iatom_o = head_array[ich]
            iatom_l = ichatbb[iatom_o]
            while True:
                ilength = 1
                while iatom_l >= 0:
                    r[0] = self._coords_unwrap[iatom_l, 0] - self._coords_unwrap[iatom_o, 0]
                    r[1] = self._coords_unwrap[iatom_l, 1] - self._coords_unwrap[iatom_o, 1]
                    r[2] = self._coords_unwrap[iatom_l, 2] - self._coords_unwrap[iatom_o, 2]
                    rsq = r[0]*r[0] + r[1]*r[1] +r[2] * r[2]
                    self._rsq_intdist[ilength] += rsq
                    self._rsqcount_intdist[ilength] += 1
                    if ilength == 1:
                        self._rsqlist.append(rsq)

                        with open("bbb.txt",'a') as f:
                            line = "{0:d} {1:f} {2:f} {3:d} {4:d} {5:f} {6:d}\n".format(ich, rsq, np.sqrt(rsq), iatom_l+1, iatom_o+1, self._rsq_intdist[ilength],  self._rsqcount_intdist[ilength])
                            f.write(line)
                    iatom_l = ichatbb[iatom_l]
                    ilength += 1
                if ichatbb[iatom_o]<=0:
                    break
                iatom_o = ichatbb[iatom_o]
                iatom_l = ichatbb[iatom_o]

   # #########################################################################
    def _write_isinternalchaindist(self, filename):

        with open(filename, 'w') as f:
            f.writelines("# of segments(n)      <R2>(A^2)/n     <R2>(A^2)/nl2      <R2>(A^2)    Count\n")
            f.writelines("# =========================================================================\n")

            for isegments in range(1, self._nbonds_bb_max+1):
                self._rsqavg_intdist[isegments] = self._rsq_intdist[isegments] /self._rsqcount_intdist[isegments]
                lines = "{0:8d}     {1:12.3f}         {2:12.3f}     {3:12.3f}      {4:10d}\n".format(isegments,
                                                                                   self._rsqavg_intdist[isegments]/isegments,
                                                                                   self._rsqavg_intdist[isegments] / isegments/ self._lavg2,
                                                                                   self._rsqavg_intdist[isegments],
                                                                                   self._rsqcount_intdist[isegments])
                f.writelines(lines)

    # #########################################################################
    # def odf_intra_cython(self, initframe=0, endframe=None, stride=1, filename="odf_intra.dat"):
    #
    #     """
    #     Calculate intra chain orientation correlation (for persistence length)
    #
    #     Calculate orientational correlation functions
    #     for intra-chain vectors as function of chemical distance.
    #     It also calculates components and 4th moments.
    #     The input is supposed to contain a multiple of nvec vectors.
    #
    #     """
    #
    #     nmols_array, neigbors = self._trajectory.topology.get_array_mols_neigh()
    #     nchains = len(nmols_array)
    #     all_bb_bonds, nbonds_bb_perch = self._trajectory.topology.get_all_bb_bonds()
    #     try:
    #         nbonds_bb_max = max(nbonds_bb_perch.values())
    #     except:
    #         nbonds_bb_max = 1
    #     coords_t0_wrapped = self._trajectory.universe.trajectory[0].positions
    #     box_dimensions = self._trajectory.universe.trajectory[0].dimensions[0:3]
    #     iatch = self._trajectory.topology._iatch
    #     nframes = self._trajectory.get_numframes()
    #     if endframe is None:
    #         endframe = nframes
    #
    #     pag.setup_odf_intra(nbonds_bb_max)
    #     uux = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
    #     uuy = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
    #     uuz = np.zeros((nchains, nbonds_bb_max), dtype=np.float32)
    #     for iframe in range(initframe ,endframe, stride):
    #         if iframe %100 == 0:
    #             print("{} of {}".format(iframe, nframes ))
    #         coords_t0_wrapped = self._trajectory.universe.trajectory[iframe].positions
    #         box_dimensions = self._trajectory.universe.trajectory[iframe].dimensions[0:3]
    #         coords_unwrap = pag.unwrap(coords_t0_wrapped, nmols_array, neigbors,
    #                                    box_dimensions, iframe=iframe)
    #
    #         pag.odf_intra(iframe, nchains, all_bb_bonds, coords_unwrap, iatch, uux, uuy, uuz)
    #
    #     pag.avg_write_odf_intra(nframes, nbonds_bb_max, filename)






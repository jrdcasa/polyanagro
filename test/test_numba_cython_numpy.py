import unittest
import datetime
import utils
import topology


# def worker(trj, q_in, q_out):
#
#
#     input = q_in.get()
#
#     coords_t0_wrapped = input[0]
#     nmols_array = input[1]
#     l_neigh_array = input[2]
#     box_dimensions = input[3]
#     iframe = input[4]
#

#     q_out.put(coords_unwrap_numba)

def worker(iproc, coords_t0_wrapped, box_dimensions, nmols_array, l_neigh_array, q_out=None):

    nblocks =len(coords_t0_wrapped)

    for iblock in range(nblocks):
        iframe = coords_t0_wrapped[iblock][0]
        c = coords_t0_wrapped[iblock][1]
        b = box_dimensions[iblock][1]
        coords_unwrap_numba = top.unwrap_numba(c, nmols_array, l_neigh_array,
                                               b, iframe=iframe)
        if q_out:
            q_out.put(coords_unwrap_numba[0:500,:])


class test_NumbaCythonNumpy(unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        cls.log = utils.init_logger("Output", fileoutput="NumbaCythonNumpy/test_performance.log", append=False, inscreen=False)
        cls.log.info("\n\tJob Information\n\t---------------")
        m = "\n\t***************** START Numba-Cython-Python TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\tStarting: \t {}\n".format(now))

    # #########################################################################
    def setup_trj(self):

        xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
        xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
        xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
        filename_tpr = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
        filename_psf = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/00-GENERATE/namd_out.psf"
        self.stride = 1

        # Setup Trajectory
        self.trj = topology.ExtTrajectory([xtc1, xtc2, xtc3], logger=self.log)
        #trj_psf = topology.Trajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

        # Setup topology
        top = topology.Topology(logger=self.log)
        top.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
        self.trj.set_top_universe(top)

    # #########################################################################
    def setup_trj_small(self):

        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
        #filename_psf = "../data/0003Ch-C020-002br04/namd_out.psf"

        self.stride = 1

        # Setup Trajectory
        self.trj = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_tpr, logger=self.log)

        #self.trj_small_psf = top.Trajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

    # #########################################################################
    def setup_trj_hd(self):

        path = "/media/jramos/PEBranch/0100Ch-C500-009br02/T-450K/06-NPT-100ns-Parrinello/"
        xtc1 = path + "01-RESTART-0000-0500ns/traj_comp.xtc"
        xtc2 = path + "02-RESTART-0500-1000ns/traj_comp.total.xtc"
        filename_tpr =  path + "01-RESTART-0000-0500ns/topol.tpr"

        self.stride = 1

        # Setup Trajectory
        self.trj = topology.ExtTrajectory([xtc1, xtc2], topfile=filename_tpr, logger=self.log)
        #self.trj_small_psf = topology.Trajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

     # #########################################################################
    # def test_01_unwrap_rg(self):
    #
    #     """
    #     	@@@@@@@@ Time setup trajectory: 	 28.26 seconds
    #     	@@@@@@@@ Time Unwrap coordinates Cython: 	 8.02 seconds
    #         @@@@@@@@ Time Unwrap coordinates PurePython: 	 73.78 seconds
    #         @@@@@@@@ Time Unwrap coordinates Numba: 	 8.94 seconds
    #         @@@@@@@@ Time Unwrap coordinates Numba nopython: 	 180.51 seconds
    #
    #         @@@@@@@@ Total time Test01: 	 299.52 seconds
    #     """
    #
    #     # START TIME **********
    #     start_t01 = datetime.datetime.now()
    #     m = "\n\t================================ START TEST_01 ================================\n"+ \
    #         "                 Unwrapping coordinates, Rg calculation  \n"
    #     print(m) if self.log is None else self.log.info(m)
    #
    #     # SETUP TRAJECTORY **********
    #     start_t01a = datetime.datetime.now()
    #
    #     #self.setup_trj()
    #     self.setup_trj_small()
    #     #self.setup_trj_hd()
    #
    #     e1 = datetime.datetime.now()-start_t01a
    #     m =  "\n\t@@@@@@@@ Time setup trajectory: \t {0:.2f} seconds\n".format(float(e1.total_seconds()))
    #     print(m) if self.log is None else self.log.info(m)
    #
    #     ini = 0
    #     nframes = self.trj.nframes
    #     stride = 1
    #     isrg = True
    #     # UNWRAP COORDINATES, Rg  CYTHON**********
    #     start_t01b = datetime.datetime.now()
    #     start_t01d = datetime.datetime.now()
    #     freq_out = (nframes - 1) * 0.05
    #     if freq_out <= 0:
    #        freq_out = 1
    #     m = "\n\t@@@@@@@@ Cython @@@@@@@@\n"
    #     print(m) if self.log is None else self.log.info(m)
    #     mass = np.array(self.trj.topology._mass, dtype=np.float32)
    #     nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
    #     nchains = len(nmols_array)
    #     rgsq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)
    #     rg2_frame = {}
    #     for iframe in range(ini, nframes, stride):
    #
    #         if iframe%freq_out == 0:
    #             elapsed_time = datetime.datetime.now() - start_t01d
    #             m = "\tIFRAME: {0:d} in {1:s} seconds".format \
    #                 (iframe, str(elapsed_time.total_seconds()))
    #             print(m) if self.log  is None else self.log .info(m)
    #
    #         coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
    #         box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
    #         coords_unwrap_cython = top.unwrap(coords_t0_wrapped, nmols_array, l_neigh_array,
    #                                           box_dimensions, iframe=iframe)
    #         if isrg:
    #             top.calc_rg_openmp(nmols_array, mass, coords_unwrap_cython, rgsq_ich_iframe)
    #             # Calculate the average
    #             rg2_avg = 0.0
    #             rg2_std = 0.0
    #
    #             # Calculate the avaerage
    #             rg2_avg = np.mean(rgsq_ich_iframe[:, 3])
    #             # Calculate the standard deviation
    #             rg2_std = np.std(rgsq_ich_iframe[:, 3])
    #             rg2_frame[iframe] = [rg2_avg, rg2_std]
    #
    #     e2 = datetime.datetime.now()-start_t01b
    #     m =  "\n\t@@@@@@@@ Time Unwrap coordinates Cython: \t {0:.2f} seconds\n".format(float(e2.total_seconds()))
    #     print(m) if self.log is None else self.log.info(m)
    #
    #
    #     # UNWRAP COORDINATES PURE PYTHON**********
    #     start_t01c = datetime.datetime.now()
    #     nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
    #     for iframe in range(ini, nframes, stride):
    #         coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
    #         box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
    #         coords_unwrap_purepython = top.unwrap_purepython(coords_t0_wrapped, self.trj.topology,
    #                                                          box_dimensions, iframe=iframe)
    #     e3 = datetime.datetime.now()-start_t01c
    #     m =  "\t@@@@@@@@ Time Unwrap coordinates PurePython: \t {0:.2f} seconds\n".format(float(e3.total_seconds()))
    #     print(m) if self.log is None else self.log.info(m)
    #
    #     # UNWRAP COORDINATES NUMBA**********
    #     start_t01d = datetime.datetime.now()
    #     freq_out = (nframes-1) * 0.05
    #     m = "\n\t@@@@@@@@ Numba @@@@@@@@\n"
    #     print(m) if self.log is None else self.log.info(m)
    #     mass = np.array(self.trj.topology._mass, dtype=np.float32)
    #     nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
    #     nchains = len(nmols_array)
    #     rgsq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)
    #     rg2_frame_numba = {}
    #     for iframe in range(ini, nframes, stride):
    #
    #         if iframe%freq_out == 0:
    #             elapsed_time = datetime.datetime.now() - start_t01d
    #             m = "\tIFRAME: {0:d} in {1:s} seconds".format \
    #                 (iframe, str(elapsed_time.total_seconds()))
    #             print(m) if self.log  is None else self.log .info(m)
    #
    #         coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
    #         box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
    #         coords_unwrap_numba = top.unwrap_numba(coords_t0_wrapped, nmols_array, l_neigh_array,
    #                                                          box_dimensions, iframe=iframe)
    #         if isrg:
    #             top.calc_rg_openmp(nmols_array, mass, coords_unwrap_numba, rgsq_ich_iframe)
    #             # Calculate the average
    #             rg2_avg = 0.0
    #             rg2_std = 0.0
    #
    #             # Calculate the avaerage
    #             rg2_avg = np.mean(rgsq_ich_iframe[:, 3])
    #             # Calculate the standard deviation
    #             rg2_std = np.std(rgsq_ich_iframe[:, 3])
    #             rg2_frame_numba[iframe] = [rg2_avg, rg2_std]
    #
    #     e4 = datetime.datetime.now()-start_t01d
    #     m =  "\t@@@@@@@@ Time Unwrap coordinates Numba: \t {0:.2f} seconds\n".format(float(e4.total_seconds()))
    #     print(m) if self.log is None else self.log.info(m)
    #
    #     # UNWRAP COORDINATES NUMBA PYTHON FUNC**********
    #     start_t01d = datetime.datetime.now()
    #     for iframe in range(ini, nframes, stride):
    #         coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
    #         box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
    #         nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
    #         coords_unwrap_numba_func = top.unwrap_numba(coords_t0_wrapped, nmols_array, l_neigh_array,
    #                                                          box_dimensions, iframe=iframe, numba=False)
    #     e5 = datetime.datetime.now()-start_t01d
    #     m =  "\t@@@@@@@@ Time Unwrap coordinates Numba nopython: \t {0:.2f} seconds\n".format(float(e5.total_seconds()))
    #     print(m) if self.log is None else self.log.info(m)
    #
    #     # Compare arrays
    #     natoms = coords_t0_wrapped.shape[0]
    #     print("Iframe: ", iframe)
    #     for i in range(natoms):
    #         for j in range(3):
    #             if (coords_unwrap_cython[i,j] - coords_unwrap_numba[i,j]) > 1e-3:
    #                 print(i, j, coords_unwrap_cython[i,j], coords_unwrap_numba[i,j], box_dimensions )
    #
    #
    #     # np.testing.assert_array_almost_equal(coords_unwrap_cython, coords_unwrap_purepython, decimal=3)
    #     np.testing.assert_array_almost_equal(coords_unwrap_cython, coords_unwrap_numba, decimal=3)
    #     # np.testing.assert_array_almost_equal(coords_unwrap_cython, coords_unwrap_numba_func, decimal=3)
    #
    #     # END TIME **********
    #     end_t01 = datetime.datetime.now()
    #     e = end_t01 -start_t01
    #     m =  "\n\t@@@@@@@@ Total time Test01: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
    #     m += "\t============== END   TEST_01 ================================"
    #     print(m) if self.log is None else self.log.info(m)
    #
    # # #########################################################################
    # # def test_02_parallel_by_frames(self):
    # #
    # #     # START TIME **********
    # #     start_t01 = datetime.datetime.now()
    # #     m = "\n\t================================ START TEST_02 ================================\n"+ \
    # #         "                 Multiprocessing by frame  \n"
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     # SETUP TRAJECTORY **********
    # #     start_t01a = datetime.datetime.now()
    # #
    # #     #self.setup_trj_small()
    # #     self.setup_trj()
    # #     #self.setup_trj_hd()
    # #
    # #     e1 = datetime.datetime.now()-start_t01a
    # #     m =  "\n\t@@@@@@@@ Time setup trajectory: \t {0:.2f} seconds\n".format(float(e1.total_seconds()))
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     ini = 0
    # #     nframes = self.trj.nframes
    # #     stride = 1
    # #
    # #     # UNWRAP COORDINATES NUMBA**********
    # #     start_t01d = datetime.datetime.now()
    # #     freq_out = (nframes-1) * 0.05
    # #      if freq_out <= 0:
    # #        freq_out = 1
    # #     m = "\n\t@@@@@@@@ Numba @@@@@@@@\n"
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     ncpus = mp.cpu_count()
    # #     print(f'starting computations on {ncpus} cores')
    # #
    # #     iframe = ini
    # #     end = nframes
    # #     block_size = 10
    # #     if (nframes - ini) < block_size:
    # #        block_size = nframes - ini
    # #     nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
    # #     ncpus = 6
    # #     while iframe < end:
    # #
    # #         if iframe%freq_out == 0:
    # #             elapsed_time = datetime.datetime.now() - start_t01d
    # #             m = "\tIFRAME: {0:d} in {1:s} seconds".format \
    # #                  (iframe, str(elapsed_time.total_seconds()))
    # #             print(m) if self.log  is None else self.log .info(m)
    # #
    # #         processes = []
    # #         q_out = mp.Queue()
    # #         for iproc in range(ncpus):
    # #             try:
    # #                 coords_t0_wrapped = []
    # #                 box_dimensions = []
    # #                 for iblock in range(block_size):
    # #                     c = self.trj.universe.trajectory[iframe].positions
    # #                     coords_t0_wrapped.append([iframe,c])
    # #                     b = self.trj.universe.trajectory[iframe].dimensions[0:3]
    # #                     box_dimensions.append([iframe, b])
    # #                     iframe += 1
    # #                 p = mp.Process(target=worker, args=(iframe, iproc, coords_t0_wrapped, box_dimensions, nmols_array, l_neigh_array, q_out))
    # #                 processes.append(p)
    # #             except IndexError:
    # #                 p = mp.Process(target=worker, args=(iframe, iproc, coords_t0_wrapped, box_dimensions, nmols_array, l_neigh_array, q_out))
    # #                 processes.append(p)
    # #                 break
    # #
    # #             iframe += stride - 1
    # #             if iframe >= end: break
    # #
    # #         for p in processes:
    # #             p.start()
    # #
    # #         for p in processes:
    # #             p.join()
    # #         print(q_out.get())
    # #
    # #     e4 = datetime.datetime.now()-start_t01d
    # #     m =  "\t@@@@@@@@ Time Unwrap coordinates Numba: \t {0:.2f} seconds\n".format(float(e4.total_seconds()))
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     # END TIME **********
    # #     end_t01 = datetime.datetime.now()
    # #     e = end_t01 -start_t01
    # #     m =  "\n\t@@@@@@@@ Total time Test01: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
    # #     m += "\t================================ END   TEST_02 ================================"
    # #     print(m) if self.log is None else self.log.info(m)
    # #     now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    # #     self.log.info("\t\tEnd: \t {}\n".format(now))
    # #     print("H")
    #
    # # #########################################################################
    # # def test_03_same_results_serialnumba_parallel(self):
    # #
    # #     # START TIME **********
    # #     start_t01 = datetime.datetime.now()
    # #     m = "\n\t================================ START TEST_03 ================================\n"+ \
    # #         "             Multiprocessing and numba in serial give the same results  \n"
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     # SETUP TRAJECTORY **********
    # #     start_t01a = datetime.datetime.now()
    # #
    # #     # self.setup_trj_small()
    # #     # self.setup_trj()
    # #     self.setup_trj_hd()
    # #
    # #     e1 = datetime.datetime.now()-start_t01a
    # #     m =  "\n\t@@@@@@@@ Time setup trajectory: \t {0:.2f} seconds\n".format(float(e1.total_seconds()))
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     ini = 100
    # #     nframes = 101
    # #     stride = 1
    # #
    # #     # UNWRAP COORDINATES NUMBA IN SERIAL**********
    # #     start_t01d = datetime.datetime.now()
    # #     freq_out = (nframes-1) * 0.05
    # #     m = "\n\t@@@@@@@@ Numba @@@@@@@@\n"
    # #     print(m) if self.log is None else self.log.info(m)
    # #     for iframe in range(ini, nframes, stride):
    # #
    # #         if iframe%freq_out == 0:
    # #             elapsed_time = datetime.datetime.now() - start_t01d
    # #             m = "\tIFRAME: {0:d} in {1:s} seconds".format \
    # #                 (iframe, str(elapsed_time.total_seconds()))
    # #             print(m) if self.log  is None else self.log .info(m)
    # #
    # #         coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
    # #         box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
    # #         nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
    # #         coords_unwrap_numba = top.unwrap_numba(coords_t0_wrapped, nmols_array, l_neigh_array,
    # #                                                box_dimensions, iframe=iframe)
    # #
    # #     e4 = datetime.datetime.now()-start_t01d
    # #     m =  "\t@@@@@@@@ Time Unwrap coordinates Numba: \t {0:.2f} seconds\n".format(float(e4.total_seconds()))
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     # UNWRAP COORDINATES NUMBA IN PARALLEL**********
    # #     start_t01d = datetime.datetime.now()
    # #     freq_out = (nframes-1) * 0.05
    # #     if freq_out <= 0:
    # #         freq_out = 1
    # #     m = "\n\t@@@@@@@@ Numba @@@@@@@@\n"
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     ncpus = mp.cpu_count()
    # #     print(f'starting computations on {ncpus} cores')
    # #
    # #     ini = 100
    # #     nframes = 101
    # #     stride = 1
    # #     iframe = ini
    # #     end = nframes
    # #     block_size = 10
    # #     if (nframes - ini) < block_size:
    # #         block_size = nframes - ini
    # #     nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
    # #     ncpus = 6
    # #     while iframe < end:
    # #
    # #         if iframe%freq_out == 0:
    # #             elapsed_time = datetime.datetime.now() - start_t01d
    # #             m = "\tIFRAME: {0:d} in {1:s} seconds".format \
    # #                  (iframe, str(elapsed_time.total_seconds()))
    # #             print(m) if self.log  is None else self.log .info(m)
    # #
    # #         processes = []
    # #         q_out = mp.Queue()
    # #         for iproc in range(ncpus):
    # #             try:
    # #                 coords_t0_wrapped = []
    # #                 box_dimensions = []
    # #                 for iblock in range(block_size):
    # #                     c = self.trj.universe.trajectory[iframe].positions
    # #                     coords_t0_wrapped.append([iframe,c])
    # #                     b = self.trj.universe.trajectory[iframe].dimensions[0:3]
    # #                     box_dimensions.append([iframe, b])
    # #                     iframe += 1
    # #                 p = mp.Process(target=worker, args=(iproc, coords_t0_wrapped, box_dimensions, nmols_array, l_neigh_array, q_out))
    # #                 processes.append(p)
    # #             except IndexError:
    # #                 p = mp.Process(target=worker, args=(iframe, iproc, coords_t0_wrapped, box_dimensions, nmols_array, l_neigh_array, q_out))
    # #                 processes.append(p)
    # #                 break
    # #
    # #             iframe += stride - 1
    # #             if iframe >= end: break
    # #
    # #         for p in processes:
    # #             p.start()
    # #
    # #         for p in processes:
    # #             p.join()
    # #
    # #         c_parallel = q_out.get()
    # #
    # #     e4 = datetime.datetime.now()-start_t01d
    # #     m =  "\t@@@@@@@@ Time Unwrap coordinates Numba: \t {0:.2f} seconds\n".format(float(e4.total_seconds()))
    # #     print(m) if self.log is None else self.log.info(m)
    # #
    # #     np.testing.assert_array_almost_equal(coords_unwrap_numba[0:500,:], c_parallel, decimal=3)
    # #
    # #     # END TIME **********
    # #     end_t01 = datetime.datetime.now()
    # #     e = end_t01 -start_t01
    # #     m =  "\n\t@@@@@@@@ Total time Test01: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
    # #     m += "\t================================ END   TEST_02 ================================"
    # #     print(m) if self.log is None else self.log.info(m)
    # #     now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
    # #     self.log.info("\t\tEnd: \t {}\n".format(now))
    # #     print("H")

    # #########################################################################
    def test_03_calc_distances(self):

        """
        	@@@@@@@@ Time setup trajectory: 	 28.26 seconds
        	@@@@@@@@ Time Unwrap coordinates Cython: 	 8.02 seconds
            @@@@@@@@ Time Unwrap coordinates PurePython: 	 73.78 seconds
            @@@@@@@@ Time Unwrap coordinates Numba: 	 8.94 seconds
            @@@@@@@@ Time Unwrap coordinates Numba nopython: 	 180.51 seconds

            @@@@@@@@ Total time Test01: 	 299.52 seconds
        """

        # START TIME **********
        start_t01 = datetime.datetime.now()
        m = "\n\t================================ START TEST_03 ================================\n"+ \
            "                 Calculate distances from serial and openmp  \n"
        print(m) if self.log is None else self.log.info(m)

        # SETUP TRAJECTORY **********
        start_t01a = datetime.datetime.now()

        self.setup_trj()
        #self.setup_trj_small()
        #self.setup_trj_hd()

        e1 = datetime.datetime.now()-start_t01a
        # m =  "\n\t@@@@@@@@ Time setup trajectory: \t {0:.2f} seconds\n".format(float(e1.total_seconds()))
        # print(m) if self.log is None else self.log.info(m)
        #
        ini = 0
        nframes = self.trj.nframes
        stride = 1

        # # UNWRAP COORDINATES, Rg  CYTHON**********
        # start_t01b = datetime.datetime.now()
        # start_t01d = datetime.datetime.now()
        # freq_out = (nframes - 1) * 0.05
        # if freq_out <= 0:
        #    freq_out = 1
        # m = "\n\t@@@@@@@@ Cython @@@@@@@@\n"
        # print(m) if self.log is None else self.log.info(m)
        # mass = np.array(self.trj.topology._mass, dtype=np.float32)
        # nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
        # nchains = len(nmols_array)
        # rgsq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)
        # rg2_frame = {}
        # for iframe in range(ini, nframes, stride):
        #
        #     if iframe%freq_out == 0:
        #         elapsed_time = datetime.datetime.now() - start_t01d
        #         m = "\tIFRAME: {0:d} in {1:s} seconds".format \
        #             (iframe, str(elapsed_time.total_seconds()))
        #         print(m) if self.log  is None else self.log .info(m)
        #
        #     coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions


        #     box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
        #     coords_unwrap_cython = top.unwrap(coords_t0_wrapped, nmols_array, l_neigh_array,
        #                                       box_dimensions, iframe=iframe)
        #     if isrg:
        #         top.calc_rg_openmp(nmols_array, mass, coords_unwrap_cython, rgsq_ich_iframe)
        #         # Calculate the average
        #         rg2_avg = 0.0
        #         rg2_std = 0.0
        #
        #         # Calculate the avaerage
        #         rg2_avg = np.mean(rgsq_ich_iframe[:, 3])
        #         # Calculate the standard deviation
        #         rg2_std = np.std(rgsq_ich_iframe[:, 3])
        #         rg2_frame[iframe] = [rg2_avg, rg2_std]
        #
        # e2 = datetime.datetime.now()-start_t01b
        # m =  "\n\t@@@@@@@@ Time Unwrap coordinates Cython: \t {0:.2f} seconds\n".format(float(e2.total_seconds()))
        # print(m) if self.log is None else self.log.info(m)
        #
        #
        # # UNWRAP COORDINATES PURE PYTHON**********
        # start_t01c = datetime.datetime.now()
        # nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
        # for iframe in range(ini, nframes, stride):
        #     coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
        #     box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
        #     coords_unwrap_purepython = top.unwrap_purepython(coords_t0_wrapped, self.trj.topology,
        #                                                      box_dimensions, iframe=iframe)
        # e3 = datetime.datetime.now()-start_t01c
        # m =  "\t@@@@@@@@ Time Unwrap coordinates PurePython: \t {0:.2f} seconds\n".format(float(e3.total_seconds()))
        # print(m) if self.log is None else self.log.info(m)
        #
        # # UNWRAP COORDINATES NUMBA**********
        # start_t01d = datetime.datetime.now()
        # freq_out = (nframes-1) * 0.05
        # m = "\n\t@@@@@@@@ Numba @@@@@@@@\n"
        # print(m) if self.log is None else self.log.info(m)
        # mass = np.array(self.trj.topology._mass, dtype=np.float32)
        # nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
        # nchains = len(nmols_array)
        # rgsq_ich_iframe = np.zeros([nchains, 4], dtype=np.float32)
        # rg2_frame_numba = {}
        # for iframe in range(ini, nframes, stride):
        #
        #     if iframe%freq_out == 0:
        #         elapsed_time = datetime.datetime.now() - start_t01d
        #         m = "\tIFRAME: {0:d} in {1:s} seconds".format \
        #             (iframe, str(elapsed_time.total_seconds()))
        #         print(m) if self.log  is None else self.log .info(m)
        #
        #     coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
        #     box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
        #     coords_unwrap_numba = top.unwrap_numba(coords_t0_wrapped, nmols_array, l_neigh_array,
        #                                                      box_dimensions, iframe=iframe)
        #     if isrg:
        #         top.calc_rg_openmp(nmols_array, mass, coords_unwrap_numba, rgsq_ich_iframe)
        #         # Calculate the average
        #         rg2_avg = 0.0
        #         rg2_std = 0.0
        #
        #         # Calculate the avaerage
        #         rg2_avg = np.mean(rgsq_ich_iframe[:, 3])
        #         # Calculate the standard deviation
        #         rg2_std = np.std(rgsq_ich_iframe[:, 3])
        #         rg2_frame_numba[iframe] = [rg2_avg, rg2_std]
        #
        # e4 = datetime.datetime.now()-start_t01d
        # m =  "\t@@@@@@@@ Time Unwrap coordinates Numba: \t {0:.2f} seconds\n".format(float(e4.total_seconds()))
        # print(m) if self.log is None else self.log.info(m)
        #
        # # UNWRAP COORDINATES NUMBA PYTHON FUNC**********
        # start_t01d = datetime.datetime.now()
        # for iframe in range(ini, nframes, stride):
        #     coords_t0_wrapped = self.trj.universe.trajectory[iframe].positions
        #     box_dimensions = self.trj.universe.trajectory[iframe].dimensions[0:3]
        #     nmols_array, l_neigh_array = self.trj.topology.get_array_mols_neigh()
        #     coords_unwrap_numba_func = top.unwrap_numba(coords_t0_wrapped, nmols_array, l_neigh_array,
        #                                                      box_dimensions, iframe=iframe, numba=False)
        # e5 = datetime.datetime.now()-start_t01d
        # m =  "\t@@@@@@@@ Time Unwrap coordinates Numba nopython: \t {0:.2f} seconds\n".format(float(e5.total_seconds()))
        # print(m) if self.log is None else self.log.info(m)
        #
        # # Compare arrays
        # natoms = coords_t0_wrapped.shape[0]
        # print("Iframe: ", iframe)
        # for i in range(natoms):
        #     for j in range(3):
        #         if (coords_unwrap_cython[i,j] - coords_unwrap_numba[i,j]) > 1e-3:
        #             print(i, j, coords_unwrap_cython[i,j], coords_unwrap_numba[i,j], box_dimensions )
        #
        #
        # # np.testing.assert_array_almost_equal(coords_unwrap_cython, coords_unwrap_purepython, decimal=3)
        # np.testing.assert_array_almost_equal(coords_unwrap_cython, coords_unwrap_numba, decimal=3)
        # # np.testing.assert_array_almost_equal(coords_unwrap_cython, coords_unwrap_numba_func, decimal=3)

        # END TIME **********
        end_t01 = datetime.datetime.now()
        e = end_t01 - start_t01
        m =  "\n\t@@@@@@@@ Total time Test03: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)
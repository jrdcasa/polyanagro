import unittest
import datetime
import numpy as np
import polyanagro as pag
import gc
import multiprocessing

class InternalCoordinatesTests (unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        cls.log = pag.init_logger("Output", fileoutput="IC_test/test_internalcoordinates.log", append=False, inscreen=False)
        cls.log.info("\n\tJob Information\n\t---------------")
        m = "\n\t\t***************** START Internal Coordinates TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

        # Create the trajectory and topology to test the internal coordinates.
        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
        cls.trj01 = pag.ExtTrajectory([xtc1], topfile=filename_tpr, logger=cls.log)

    # #########################################################################
    def aux_distances(self, ref, conf, itera=100):

        m = "\n\t ******* DISTANCE_ARRAY NUMPY *******\n"
        start_time = datetime.datetime.now()
        for _ in range(itera):
            dij_np = pag.distance_array_numpypython(ref, conf)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m += "\t Distance_array numpypython ({0:d} times) --> time: {1:s} seconds".format(itera, str(elapsed_time.total_seconds()))
        self.log.info(m)

        m = "\t ******* DISTANCE_ARRAY PUREPYTHON *******\n"
        start_time = datetime.datetime.now()
        for _ in range(itera):
            dij_pp = pag.distance_array_purepython(ref, conf)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m += "\t Distance_array purepython ({0:d} times) --> time: {1:s} seconds".format(itera, str(elapsed_time.total_seconds()))
        self.log.info(m)
        m = "\t ******* DISTANCE_ARRAY CYTHON *******\n"
        start_time = datetime.datetime.now()
        for _ in range(itera):
            dij_c = pag.distance_array(ref, conf)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m += "\t Distance_array cython ({0:d} times) --> time: {1:s} seconds".format(itera, str(elapsed_time.total_seconds()))
        self.log.info(m)

        # Compare arrays
        np.testing.assert_array_almost_equal(dij_np, dij_pp, decimal=3)
        np.testing.assert_array_almost_equal(dij_pp, dij_c,  decimal=3)

    # #########################################################################
    def aux_diagonal_distances(self, ref, conf, itera=100):

        m = "\t ******* DISTANCE_ARRAY CYTHON *******\n"
        start_time = datetime.datetime.now()
        for _ in range(itera):
            dij_c = pag.distance_diagonal_array(ref, conf)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m += "\t Distance_diagonal_array cython ({0:d} times) --> time: {1:s} seconds".format(itera, str(elapsed_time.total_seconds()))
        self.log.info(m)
        return dij_c

    # #########################################################################
    def test_01(self):

        m = "\n\t============== START TEST_01 ================================\n"+ \
            "              Bend, Improper and dihedral angles\n"
        print(m) if self.log is None else self.log.info(m)

        m = "\t3 Chains 20 backbone C and 2 C2 branches for each chain\n"
        m+= "\tTotal atoms 60 + 4*3 = 72 atoms\n"
        print(m) if self.log is None else self.log.info(m)

        m = "\tIn this test, the coordinates from Frame 0 are used.\n"
        print(m) if self.log is None else self.log.info(m)
        coords_t0_wrapped = self.trj01.universe.atoms.positions
        box_dimensions = self.trj01.universe.dimensions[0:3]

        # Unwrap coordinates Pure python
        start_time = datetime.datetime.now()
        for _ in range(5000):
            coords_t0_unwrapped01 = pag.unwrap_purepython(coords_t0_wrapped, self.trj01.topology,
                                                          box_dimensions, iframe=0 )
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\tUnwrapping coordinates pure python (5000 times) --> time: {0:s} seconds".format(str(elapsed_time.total_seconds()))
        self.log.info(m)

        # Unwrap coordinates Cython
        nmols_array, l_neigh_array = self.trj01.topology.get_array_mols_neigh()
        start_time = datetime.datetime.now()
        for _ in range(5000):
            coords_t0_unwrapped02 = pag.unwrap(coords_t0_wrapped, nmols_array, l_neigh_array,
                                               box_dimensions, iframe=0)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\tUnwrapping coordinates cython (5000 times) --> time: {0:s} seconds".format(str(elapsed_time.total_seconds()))
        self.log.info(m)

        # Unwrap and Nojump coordinates Cython
        n_image = np.zeros((self.trj01.topology._natoms, 3), dtype=np.int32)
        start_time = datetime.datetime.now()
        for _ in range(5000):
            coords_t0_unwrapped_nojump = pag.unwrap_nojump(coords_t0_wrapped,coords_t0_wrapped, n_image,
                                                           nmols_array, l_neigh_array, box_dimensions, iframe=0)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\tUnwrapping coordinates cython (5000 times) --> time: {0:s} seconds".format(str(elapsed_time.total_seconds()))
        self.log.info(m)

        # Compare arrays
        np.testing.assert_array_almost_equal(coords_t0_unwrapped01, coords_t0_unwrapped02, decimal=3)
        np.testing.assert_array_almost_equal(coords_t0_unwrapped01, coords_t0_unwrapped_nojump, decimal=3)

        # Calculate Bend angle
        c1 = coords_t0_unwrapped01[3,:]
        c2 = coords_t0_unwrapped01[4,:]
        c3 = coords_t0_unwrapped01[5,:]
        b1_d = pag.bend_angle_purepython(c1, c2, c3)
        b1_r = pag.bend_angle_purepython(c1, c2, c3, radians=True)
        vmd_value = 113.609 # degrees
        m = "\n\t ******* BEND ANGLE *******\n"
        m += "\tCoordinates Atom 3: {0:.3f} {1:.3f} {2:.3f}\n".format(c1[0],c1[1], c1[2] )
        m += "\tCoordinates Atom 4: {0:.3f} {1:.3f} {2:.3f}\n".format(c2[0], c2[1], c2[2])
        m += "\tCoordinates Atom 5: {0:.3f} {1:.3f} {2:.3f}\n".format(c3[0], c3[1], c3[2])
        m += "\t ** Bend calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(b1_d, b1_r)
        m += "\t ** Bend VMD calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(vmd_value, vmd_value*np.pi/180.0)
        self.log.info(m)

        # Calculate dihedral
        c1 = coords_t0_unwrapped01[3,:]
        c2 = coords_t0_unwrapped01[4,:]
        c3 = coords_t0_unwrapped01[5,:]
        c4 = coords_t0_unwrapped01[6,:]
        d1_d = pag.dihedral_angle_purepython(c1, c2, c3, c4)
        d1_r = pag.dihedral_angle_purepython(c1, c2, c3, c4, radians=True)
        vmd_value = 62.632 # degrees
        m = "\n\t ******* DIHEDRAL ANGLE *******\n"
        m += "\tCoordinates Atom 3: {0:.3f} {1:.3f} {2:.3f}\n".format(c1[0],c1[1], c1[2] )
        m += "\tCoordinates Atom 4: {0:.3f} {1:.3f} {2:.3f}\n".format(c2[0], c2[1], c2[2])
        m += "\tCoordinates Atom 5: {0:.3f} {1:.3f} {2:.3f}\n".format(c3[0], c3[1], c3[2])
        m += "\tCoordinates Atom 6: {0:.3f} {1:.3f} {2:.3f}\n".format(c4[0], c4[1], c4[2])
        m += "\t ** Dihedral calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(d1_d, d1_r)
        m += "\t ** Dihedral VMD calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(vmd_value, vmd_value*np.pi/180.0)
        self.log.info(m)

        # Calculate improper
        c1 = coords_t0_unwrapped01[9,:]
        c2 = coords_t0_unwrapped01[8,:]
        c3 = coords_t0_unwrapped01[10,:]
        c4 = coords_t0_unwrapped01[20,:]
        d1_d = pag.dihedral_angle_purepython(c1, c2, c3, c4)
        d1_r = pag.dihedral_angle_purepython(c1, c2, c3, c4, radians=True)
        vmd_value = 30.726 # degrees
        m = "\n\t ******* IMPROPER ANGLE *******\n"
        m += "\tCoordinates Atom 9: {0:.3f} {1:.3f} {2:.3f}\n".format(c1[0],c1[1], c1[2] )
        m += "\tCoordinates Atom 8: {0:.3f} {1:.3f} {2:.3f}\n".format(c2[0], c2[1], c2[2])
        m += "\tCoordinates Atom 10: {0:.3f} {1:.3f} {2:.3f}\n".format(c3[0], c3[1], c3[2])
        m += "\tCoordinates Atom 20: {0:.3f} {1:.3f} {2:.3f}\n".format(c4[0], c4[1], c4[2])
        m += "\t ** Improper calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(d1_d, d1_r)
        m += "\t ** Improper VMD calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(vmd_value, vmd_value*np.pi/180.0)
        self.log.info(m)

        # Calculate cos_angle_purepython
        c1 = coords_t0_unwrapped01[9,:]
        c2 = coords_t0_unwrapped01[8,:]
        c3 = coords_t0_unwrapped01[10,:]
        c4 = coords_t0_unwrapped01[20,:]
        d1_d = pag.cos_angle_purepython(c1, c2, c3, c4)
        d1_r = pag.cos_angle_purepython(c1, c2, c3, c4, radians=True)
        vmd_value = 42.221 # degrees
        m = "\n\t ******* COS_ANGLE ANGLE *******\n"
        m += "\tCoordinates Atom 9: {0:.3f} {1:.3f} {2:.3f}\n".format(c1[0], c1[1], c1[2])
        m += "\tCoordinates Atom 8: {0:.3f} {1:.3f} {2:.3f}\n".format(c2[0], c2[1], c2[2])
        m += "\tCoordinates Atom 10: {0:.3f} {1:.3f} {2:.3f}\n".format(c3[0], c3[1], c3[2])
        m += "\tCoordinates Atom 20: {0:.3f} {1:.3f} {2:.3f}\n".format(c4[0], c4[1], c4[2])
        m += "\t ** cos_angle calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(d1_d, d1_r)
        m += "\t ** cos_angle VMD calculated: {0:.3f} degrees ({1:.3f} rad) **\n".format(vmd_value, vmd_value*np.pi/180.0)
        self.log.info(m)


        m = "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_02(self):

        m = "\n\t============== START TEST_02 ================================\n"+ \
            "                    Distance array \n"
        print(m) if self.log is None else self.log.info(m)

        # All distances in system (72 atoms). Frame 0
        ref = np.ascontiguousarray(self.trj01.universe.atoms.positions, dtype=np.float64)
        conf = np.ascontiguousarray(self.trj01.universe.atoms.positions, dtype=np.float64)

        m = "\tReference array (72 elements) and conf array (72 elements)"
        self.log.info(m)

        self.aux_distances(ref, conf)

        m = "\t============== END   TEST_02 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_03(self):

        m = "\n\t============== START TEST_03 ================================\n"+ \
            "         Distance array. Ref and conf different sizes \n"
        print(m) if self.log is None else self.log.info(m)

        # All distances in system (72 atoms). Frame 0
        ref = np.ascontiguousarray(self.trj01.universe.atoms.positions[0:10], dtype=np.float64)
        conf = np.ascontiguousarray(self.trj01.universe.atoms.positions[40:71], dtype=np.float64)

        m = "\tReference array (10 elements) and conf array (31 elements)"
        self.log.info(m)

        self.aux_distances(ref, conf)

        m = "\t============== END   TEST_03 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_04(self):

        m = "\n\t============== START TEST_04 ================================\n"+ \
            "              Diagonal Distance array \n"
        print(m) if self.log is None else self.log.info(m)

        # All distances in system (72 atoms). Frame 0
        ref = np.ascontiguousarray(self.trj01.universe.atoms.positions[0:2], dtype=np.float64)
        conf = np.ascontiguousarray(self.trj01.universe.atoms.positions[3:5], dtype=np.float64)

        m = "\tDistances between atoms 0-3 and 1-4 "
        self.log.info(m)

        d = self.aux_diagonal_distances(ref, conf, itera=250000)

        np.testing.assert_array_almost_equal(d[0], np.array([24.080, 3.183]), decimal=3)

        m = "\t============== END   TEST_04 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_05(self):

        m = "\n\t============== START TEST_05 ================================\n"+ \
            "                  Center of geometry \n"
        print(m) if self.log is None else self.log.info(m)

        itera = 1000

        #Unwrap coordinates
        coords_t0_wrapped = self.trj01.universe.atoms.positions
        box_dimensions = self.trj01.universe.atoms.dimensions[0:3]
        nmols_array, l_neigh_array = self.trj01.topology.get_array_mols_neigh()
        coords_t0_unwrapped = pag.unwrap(coords_t0_wrapped, nmols_array, l_neigh_array,
                                               box_dimensions, iframe=0)

        # Center of geometry pure python
        nchains = len(nmols_array)
        m = "\t Number of chains: {}".format(nchains)
        print(m) if self.log is None else self.log.info(m)
        cog_pp = np.zeros((nchains,3))
        cog_vmd = np.zeros((nchains,3))
        cog_vmd[0] = [23.809 ,24.668, 25.432]
        cog_vmd[1] = [15.965, 15.226,  0.824]
        cog_vmd[2] = [16.933, 15.715,  7.612]
        start_time = datetime.datetime.now()
        for _ in range(itera):
            for ich in range(nchains):
                nat = sum(map(lambda x: x>=0, nmols_array[ich]))
                c_aux = np.zeros((nat,3))
                for iat in range(nat):
                    c_aux[iat,:] = coords_t0_unwrapped[nmols_array[ich][iat],:]
                cog_pp[ich] = pag.center_of_geom_purepython(c_aux)
        end_time = datetime.datetime.now()
        elapsed_time_py = end_time - start_time

        # Center of geometry cython
        cog_cy = np.zeros((nchains,3))
        start_time = datetime.datetime.now()
        for _ in range(itera):
            for ich in range(nchains):
                nat = sum(map(lambda x: x>=0, nmols_array[ich]))
                c_aux = np.zeros((nat,3), dtype=np.float32)
                for iat in range(nat):
                    c_aux[iat,:] = coords_t0_unwrapped[nmols_array[ich][iat],:]
                cog_cy[ich] = pag.center_of_geom(c_aux)
        end_time = datetime.datetime.now()
        elapsed_time_cy = end_time - start_time
        m = ""
        for ich in range(nchains):
            m+="\t Center of geometry of chain {0:d} (Pure Python) : {1:.3f}, {2:.3f}, {3:.3f}\n".format(ich, cog_pp[ich][0], cog_pp[ich][1], cog_pp[ich][2])
            m+="\t Center of geometry of chain {0:d} (Cython)      : {1:.3f}, {2:.3f}, {3:.3f}\n".format(ich, cog_cy[ich][0], cog_cy[ich][1], cog_cy[ich][2])
            m+="\t Center of geometry of chain {0:d} (VMD)         : {1:.3f}, {2:.3f}, {3:.3f}\n\n".format(ich, cog_vmd[ich][0],cog_vmd[ich][1], cog_vmd[ich][2])
        print(m) if self.log is None else self.log.info(m)

        m = "\t ******* TIMES *******\n"
        m += "\t Center of geometry purepython ({0:d} times) --> time: {1:s} seconds\n".format(itera, str(elapsed_time_py.total_seconds()))
        m += "\t Center of geometry     cython ({0:d} times) --> time: {1:s} seconds\n".format(itera, str(elapsed_time_cy.total_seconds()))
        print(m) if self.log is None else self.log.info(m)

        m = "\t============== END   TEST_05 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_06(self):

        m = "\n\t============== START TEST_06 ================================\n"+ \
            "                  Center of mass \n"
        print(m) if self.log is None else self.log.info(m)

        itera = 1000

        #Unwrap coordinates
        #coords_t0_wrapped = self.trj01.universe.atoms.positions
        coords_t0_wrapped = self.trj01.universe.trajectory[0].positions
        box_dimensions = self.trj01.universe.atoms.dimensions[0:3]
        nmols_array, l_neigh_array = self.trj01.topology.get_array_mols_neigh()
        coords_t0_unwrapped = pag.unwrap(coords_t0_wrapped, nmols_array, l_neigh_array,
                                               box_dimensions, iframe=0)

        # Center of mass pure python
        nchains = len(nmols_array)
        m = "\t Number of chains: {}\n".format(nchains)
        print(m) if self.log is None else self.log.info(m)
        com_pp = np.zeros((nchains,3))
        com_vmd = np.zeros((nchains,3))
        com_vmd[0] = [23.799 ,24.655, 25.427]
        com_vmd[1] = [15.988, 15.237,  0.817]
        com_vmd[2] = [16.941, 15.715,  7.626]
        start_time = datetime.datetime.now()
        for _ in range(itera):
            for ich in range(nchains):
                nat = sum(map(lambda x: x>=0, nmols_array[ich]))
                c_aux = np.zeros((nat,3))
                m_aux = np.zeros(nat)
                for iat in range(nat):
                    c_aux[iat,:] = coords_t0_unwrapped[nmols_array[ich][iat],:]
                    m_aux[iat] = self.trj01.topology._mass[nmols_array[ich][iat]]
                com_pp[ich] = pag.center_of_mass_purepython(c_aux, m_aux)
        end_time = datetime.datetime.now()
        elapsed_time_py = end_time - start_time

        # Center of mass cython
        com_cy = np.zeros((nchains,3))
        start_time = datetime.datetime.now()
        for _ in range(itera):
            for ich in range(nchains):
                nat = sum(map(lambda x: x>=0, nmols_array[ich]))
                c_aux = np.zeros((nat,3))
                m_aux = np.zeros(nat)
                for iat in range(nat):
                    c_aux[iat,:] = coords_t0_unwrapped[nmols_array[ich][iat],:]
                    m_aux[iat] = self.trj01.topology._mass[nmols_array[ich][iat]]

                com_cy[ich] = pag.center_of_mass(c_aux, m_aux)
        end_time = datetime.datetime.now()
        elapsed_time_cy = end_time - start_time

        m = ""
        for ich in range(nchains):
            m+="\t Center of mass of chain {0:d} (Pure Python) : {1:.3f}, {2:.3f}, {3:.3f}\n".format(ich, com_pp[ich][0], com_pp[ich][1], com_pp[ich][2])
            m+="\t Center of mass of chain {0:d} (Cython)      : {1:.3f}, {2:.3f}, {3:.3f}\n".format(ich, com_cy[ich][0], com_cy[ich][1], com_cy[ich][2])
            m+="\t Center of mass of chain {0:d} (VMD)         : {1:.3f}, {2:.3f}, {3:.3f}\n\n".format(ich, com_vmd[ich][0],com_vmd[ich][1], com_vmd[ich][2])
        print(m) if self.log is None else self.log.info(m)

        m = "\t ******* TIMES *******\n"
        m += "\t Center of mass purepython ({0:d} times) --> time: {1:s} seconds\n".format(itera, str(elapsed_time_py.total_seconds()))
        m += "\t Center of mass     cython ({0:d} times) --> time: {1:s} seconds\n".format(itera, str(elapsed_time_cy.total_seconds()))
        print(m) if self.log is None else self.log.info(m)

        m = "\t============== END   TEST_06 ================================"
        print(m) if self.log is None else self.log.info(m)

    # # #########################################################################
    def test_07_distancearray_benchmark(self):

        m = "\n\t============== START TEST_07 ================================\n"+ \
            "         Distance array. Serial, openmp, numba \n"
        print(m) if self.log is None else self.log.info(m)

        # Create the trajectory and topology to test the internal coordinates.
        xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
        xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
        xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
        filename_tpr = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
        trj01 = pag.ExtTrajectory([xtc1], topfile=filename_tpr, logger=self.log)

        # Frame 0
        ref = np.ascontiguousarray(trj01.universe.atoms.positions[0:5000], dtype=np.float64)
        conf = np.ascontiguousarray(trj01.universe.atoms.positions[0:5000], dtype=np.float64)

        niters = 1
        # SERIAL
        start_time = datetime.datetime.now()
        print("Running {} iterations non-parallel CYTHON-EXT version".format(niters))
        for i in range(niters):
            if i == 0:
                a1, b1, c1, d1 = pag.distance_array(ref, conf, openmp=False)
            else:
                a, b, c, d = pag.distance_array(ref, conf, openmp=False)
                # Garbage collector
                del a; del b; del c; del d
                gc.collect()
        end_time = datetime.datetime.now()
        serial_time = end_time - start_time

        #OPENMP
        start_time = datetime.datetime.now()
        print("Running {} iterations openmp CYTHON-EXT version".format(niters))
        for i in range(niters):
            if i == 0:
                a2, b2, c2, d2 = pag.distance_array(ref, conf, openmp=True)
            else:
                a, b, c, d = pag.distance_array(ref, conf, openmp=True)
                # Garbage collector
                del a; del b; del c; del d
                gc.collect()
        end_time = datetime.datetime.now()
        openmp_time = end_time - start_time

        # NUMPY
        start_time = datetime.datetime.now()
        print("Running {} iterations numpy version".format(niters))
        for i in range(niters):
            if i == 0:
                a3, b3, c3, d3 = pag.distance_array_numpypython(ref, conf)
            else:
                a, b, c, d = pag.distance_array_numpypython(ref, conf)
                # Garbage collector
                del a; del b; del c; del d
                gc.collect()
        end_time = datetime.datetime.now()
        numpy_time = end_time - start_time

        # # PURE_PYTHON
        start_time = datetime.datetime.now()
        print("Running {} iterations pure python version".format(niters))
        for i in range(niters):
            if i == 0:
                a4, b4, c4, d4 = pag.distance_array_purepython(ref, conf)
            else:
                a, b, c, d = pag.distance_array_purepython(ref, conf)
                # Garbage collector
                del a; del b; del c; del d
                gc.collect()
        end_time = datetime.datetime.now()
        purepython_time = end_time - start_time

        # NUMBA
        start_time = datetime.datetime.now()
        print("Running {} iterations numba python version".format(niters))
        for i in range(niters):
            dist = np.zeros([ref.shape[0], conf.shape[0]], dtype=np.float64)
            rijx = np.zeros([ref.shape[0], conf.shape[0]], dtype=np.float64)
            rijy = np.zeros([ref.shape[0], conf.shape[0]], dtype=np.float64)
            rijz = np.zeros([ref.shape[0], conf.shape[0]], dtype=np.float64)
            if i == 0:
                a5, b5, c5, d5 = pag.calc_distance_array_numba(ref, conf)
            else:
                a, b, c, d = pag.calc_distance_array_numba(ref, conf)
                # Garbage collector
                del a; del b; del c; del d
                gc.collect()
        end_time = datetime.datetime.now()
        numba_time = end_time - start_time


        # Testing results
        np.testing.assert_almost_equal(a1, a2, decimal=6)
        np.testing.assert_almost_equal(b1, b2, decimal=6)
        np.testing.assert_almost_equal(c1, c2, decimal=6)
        np.testing.assert_almost_equal(d1, d2, decimal=6)
        del a2; del b2; del c2; del d2
        gc.collect()
        np.testing.assert_almost_equal(a1, a3, decimal=6)
        np.testing.assert_almost_equal(b1, b3, decimal=6)
        np.testing.assert_almost_equal(c1, c3, decimal=6)
        np.testing.assert_almost_equal(d1, d3, decimal=6)
        del a3; del b3; del c3; del d3
        gc.collect()
        np.testing.assert_almost_equal(a1, a4, decimal=6)
        np.testing.assert_almost_equal(b1, b4, decimal=6)
        np.testing.assert_almost_equal(c1, c4, decimal=6)
        np.testing.assert_almost_equal(d1, d4, decimal=6)
        del a4; del b4; del c4; del d4
        gc.collect()
        np.testing.assert_almost_equal(a1, a5, decimal=6)
        np.testing.assert_almost_equal(b1, b5, decimal=6)
        np.testing.assert_almost_equal(c1, c5, decimal=6)
        np.testing.assert_almost_equal(d1, d5, decimal=6)
        del a5; del b5; del c5; del d5
        gc.collect()

        del a1; del b1; del c1; del d1
        gc.collect()
        ncpus = multiprocessing.cpu_count()

        m = "\t**** TIMES FOR {} ITERATIONS, ARRAY 5000x5000 ****\n".format(niters)
        m += "\tTime for serial   run PUREPYTHON (distance_array_purepython)   : {0:>.3f} seconds\n"\
            .format(purepython_time.total_seconds())
        m += "\tTime for serial   run NUMPY      (distance_array_numpy)        : {0:>.3f} seconds\n"\
            .format(numpy_time.total_seconds())
        m += "\tTime for serial   run CYTHON-EXT (distance_array, openmp=False): {0:>.3f} seconds\n"\
            .format(serial_time.total_seconds())
        m += "\tTime for parallel run CYTHON-EXT (distance_array, openmp=True ): {0:>.3f} seconds ({1:d} cpus)\n"\
            .format(openmp_time.total_seconds(), ncpus)
        m += "\tTime for parallel run NUMBA      (calc_distance_array_numba )  : {0:>.3f} seconds \n"\
            .format(numba_time.total_seconds(), ncpus)
        m += "\tAll implementations give the same results"
        print(m) if self.log is None else self.log.info(m)
        m = "\t============== END   TEST_07 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_08_distancediagonalarray_benchmark(self):

        m = "\n\t============== START TEST_08 ================================\n"+ \
            "         Distance diagonal array. Serial, openmp, numba \n"
        print(m) if self.log is None else self.log.info(m)

        # Create the trajectory and topology to test the internal coordinates.
        xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
        xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
        xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
        filename_tpr = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
        trj01 = pag.ExtTrajectory([xtc1], topfile=filename_tpr, logger=self.log)

        # Frame 0
        ref = np.ascontiguousarray(trj01.universe.atoms.positions[0:50000], dtype=np.float64)
        conf = np.ascontiguousarray(trj01.universe.atoms.positions[0:50000], dtype=np.float64)

        niters = 10
        # SERIAL
        start_time = datetime.datetime.now()
        for i in range(niters):
            print("Iteration {} if {}".format(i, niters))
            a, b, c, d = pag.distance_diagonal_array(ref, conf, openmp=False)
            # Garbage collector
            del a; del b; del c; del d
            gc.collect()
        end_time = datetime.datetime.now()
        serial_time = end_time - start_time
        #OPENMP
        start_time = datetime.datetime.now()
        for i in range(niters):
            print("Iteration {} if {}".format(i, niters))
            a, b, c, d = pag.distance_diagonal_array(ref, conf, openmp=True)
            # Garbage collector
            del a; del b; del c; del d
            gc.collect()
        end_time = datetime.datetime.now()
        openmp_time = end_time - start_time

        a1, b1, c1, d1 = pag.distance_diagonal_array(ref, conf, openmp=False)
        a2, b2, c2, d2 = pag.distance_diagonal_array(ref, conf, openmp=True)

        np.testing.assert_equal(a1, a2)
        np.testing.assert_equal(b1, b2)
        np.testing.assert_equal(c1, c2)
        np.testing.assert_equal(d1, d2)

        ncpus = multiprocessing.cpu_count()

        m = "\t**** TIMES FOR {} ITERATIONS, ARRAY 15000x15000 ****\n".format(niters)
        m += "\tTime for serial   run CYTHON-EXT (distance_array, openmp=False): {0:>.3f} seconds\n"\
            .format(serial_time.total_seconds())
        m += "\tTime for serial   run NUMPY      (distance_array, openmp=False): {0:>.3f} seconds\n"\
            .format(serial_time.total_seconds())
        m += "\tTime for parallel run CYTHON-EXT  (distance_array, openmp=True ): {0:>.3f} seconds ({1:d} cpus)\n"\
            .format(openmp_time.total_seconds(), ncpus)
        m += "\tBoth openmp=False and openmp=True give the same results"
        print(m) if self.log is None else self.log.info(m)

        m = "\t============== END   TEST_08 ================================"
        print(m) if self.log is None else self.log.info(m)


if __name__ == '__main__':
    unittest.main ()

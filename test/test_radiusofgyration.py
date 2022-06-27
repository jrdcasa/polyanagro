import unittest
import datetime
import io
import os
import topology
import utils
import polyanagro as pag


class RadiusofGyrationTests (unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        # Initialize log
        cls.log = utils.init_logger("Output", fileoutput="test_radiusofgyration/test_chainstatistics.log",
                                   append=False, inscreen=False)
        cls.log.info("\n\tJob Information\n\t---------------")
        m = "\n\t\t***************** START RadiusofGyration TEST *****************\n"
        m += "\t\t\tEstimated time to run the test is 2-3 minutes."
        print(m) if cls.log is None else cls.log.info(m)
        cls.start = datetime.datetime.now()
        cls.log.info("\t\tStarting: \t {}\n".format(cls.start.strftime("%d-%m-%Y %H:%M:%S")))

    # #########################################################################
    @staticmethod
    def aux_print(self, s1, s2):

        for i in range(s1._natoms):

            if not (s1._topology._orderbonds[i] == s2._topology._orderbonds[i]).all():
                print (i, s1._topology._orderbonds[i])
                print (i, s2._topology._orderbonds[i])
                print("=========")

    # #########################################################################
    def aux_smalltrj(self):

        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
        filename_psf = "../data/0003Ch-C020-002br04/namd_out.psf"

        self.stride = 1

        # Setup Trajectory
        m1 = "\t\t Trajectory XTC and TOPOLOGY TPR\n"
        m2 = "\t\t"+len(m1)*"-"
        print(m1+m2) if self.log is None else self.log.info(m1+m2)
        self.trj_small = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_tpr, logger=self.log)
        m1 = "\t\t Trajectory XTC and TOPOLOGY PSF\n"
        m2 = "\t\t"+len(m1)*"-"
        print(m1+m2) if self.log is None else self.log.info(m1+m2)

        self.trj_small_psf = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

    # #########################################################################
    def aux_comparedimensions(self, diroutput, fileref, delta=0.01):

        if diroutput[-1] != "/":
            diroutput = diroutput+"/"
        success = False

        # Rg2 ---- Check the results with reference data obtained from other program (PreMCParallel)
        f1 = io.open(diroutput+"Rg.dat")
        f2 = io.open(diroutput+fileref)
        Rg1 = []
        Rg2 = []
        for item in f1:
            if item[0] == "#": continue
            Rg1.append(float(item.split()[2]))
        for item in f2:
            if item[0] == "#": continue
            Rg2.append(float(item.split()[2]))

        Rg1Rg2 = list(zip(Rg1, Rg2))
        with open(diroutput+"Rg2_compare.dat", 'w') as f:
            f.writelines("# iframe    <Rg2(A2)>this   <Rg2(A2)>PreMC   Diff(A2)\n")
            f.writelines("=====================================================\n")
            iframe = 0
            for item in Rg1Rg2:
                a = float(item[0])
                b = float(item[1])
                l = "{0:10d}   {1:>10.2f}   {2:>10.2f}   {3:>10.2f}\n".format(iframe, a, b, a - b)
                try:
                    self.assertAlmostEqual(a, b, delta=delta)
                except AssertionError as e:
                    m = "\t{}".format(e)
                    print(m) if self.log is None else self.log.info(m)
                    f1.close()
                    f2.close()
                    return success
                f.writelines(l)
                iframe += 1

        f1.close()
        f2.close()

        # Ree2 ---- Check the results with reference data obtained other program (PreMCParallel)
        f1 = io.open(diroutput+"Ree.dat")
        f2 = io.open(diroutput+fileref)
        Ree1 = []
        Ree2 = []
        for item in f1:
            if item[0] == "#": continue
            Ree1.append(float(item.split()[2]))
        for item in f2:
            if item[0] == "#": continue
            Ree2.append(float(item.split()[3]))

        Ree1Ree2 = list(zip(Ree1, Ree2))

        with open(diroutput+"Ree2_compare.dat", 'w') as f:
            f.writelines("# iframe    <Ree2(A2)>this   <Ree2(A2)>PreMC   Diff(A2)\n")
            f.writelines("=====================================================\n")
            iframe = 0
            for item in Ree1Ree2:
                a = float(item[0])
                b = float(item[1])
                l = "{0:10d}   {1:>10.2f}   {2:>10.2f}   {3:>10.2f}\n".format(iframe, a, b, a - b)
                try:
                    self.assertAlmostEqual(a, b, delta=delta)
                except AssertionError as e:
                    m = "\t{}".format(e)
                    print(m) if self.log is None else self.log.info(m)
                    f1.close()
                    f2.close()
                    return success
                f.writelines(l)
                iframe += 1

        f1.close()
        f2.close()
        success = True
        return success

    # #########################################################################
    def aux_compareinternaldistances(self, diroutput, fileref, delta=0.01):

        if diroutput[-1] != "/":
            diroutput = diroutput+"/"
        success = False

        # Rg2 ---- Check the results with reference data obtained from other program (PreMCParallel)
        f1 = io.open(diroutput+"cn_internal_distances.dat")
        f2 = io.open(diroutput+fileref)
        Ree_n_1 = []
        Ree_n_2 = []
        for item in f1:
            if item[0] == "#": continue
            Ree_n_1.append(float(item.split()[1]))
        for item in f2:
            if item[0] == "#": continue
            Ree_n_2.append(float(item.split()[2]))

        Rg1Rg2 = list(zip(Ree_n_1, Ree_n_2))
        with open(diroutput+"Ree_internal_compare.dat", 'w') as f:
            f.writelines("# iframe    <Ree>/n(A2)>this   <Ree>/n(A2)>PreMC   Diff(A2)\n")
            f.writelines("=====================================================\n")
            iframe = 0
            for item in Rg1Rg2:
                a = float(item[0])
                b = float(item[1])
                l = "{0:10d}   {1:>10.2f}   {2:>10.2f}   {3:>10.2f}\n".format(iframe, a, b, a - b)
                try:
                    self.assertAlmostEqual(a, b, delta=delta)
                except AssertionError as e:
                    m = "\t{}".format(e)
                    print(m) if self.log is None else self.log.info(m)
                    f1.close()
                    f2.close()
                    return success
                f.writelines(l)
                iframe += 1

        f1.close()
        f2.close()

        success = True
        return success


    # # #########################################################################
    def test_01_chainstatistics_small(self):

        """
        Test 01: Small chains
        """

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_01 ================================\n"+ \
            "                  Small Trajectory Rg2 and Ree2 \n"
        print(m) if self.log is None else self.log.info(m)

        # Load small trajectory
        self.aux_smalltrj()
        natoms = self.trj_small.topology.natoms

        # Create a list with atoms to calculate the end-to-end distances
        nmols_array, neigbors = self.trj_small.topology.get_array_mols_neigh()
        nchains = len(nmols_array)
        listend2end = []
        natch = 24
        natbb = 20
        for ich in range(nchains):
            ihead = 0 + ich*natch
            itail = (natbb-1) + (ich*natch)
            listend2end.append([ich, ihead, itail])

        # Backbone atoms. In this trajectory the branches are located at the end of the chain
        backbone_list_atoms = []
        is_bb_atoms = natoms * [False]
        for ich in range(nchains):

            ihead = 0 + ich * natch
            l = []
            for ibb in range(ihead, ihead+natbb):
                l.append(ibb)
                is_bb_atoms[ibb] = True
            backbone_list_atoms.append(l)

        # Calculations Cn not calc_distances                          ----> Test 01a ==========
        start_t01a = datetime.datetime.now()
        objcalc = pag.Chain_Statistics(self.trj_small, dt=0.1,
                                       stride=self.stride, log=self.log)
        objcalc.calculate(listend2end, diroutput="test_radiusofgyration/test01/",
                          acfE2E=True, backbone_list_atoms=backbone_list_atoms,
                          isbbatom=is_bb_atoms, calc_Cn_bonds_distances=False)
        e1 = datetime.datetime.now()-start_t01a
        m =  "\tTotal time Test01a (Cn without calc_distances): \t {0:.2f} seconds\n".\
            format(float(e1.total_seconds()))
        m1 = "\t"+len(m)*":"+"\n"
        print(m1+m+m1) if self.log is None else self.log.info(m1+m+m1)
        os.rename("test_radiusofgyration/test01/Cn.dat", "test_radiusofgyration/test01/Cn_withoutdist.dat")

        # # Calculations Cn calc_distances                            ----> Test 01b ==========
        # start_t01b = datetime.datetime.now()
        # objcalc = pag.Chain_Statistics(self.trj_small, dt=0.1,
        #                                stride=self.stride, log=self.log)
        # objcalc.calculate(listend2end,diroutput="test_radiusofgyration/test01/",
        #                   acfE2E=True, backbone_list_atoms=backbone_list_atoms,
        #                   isbbatom=is_bb_atoms, calc_Cn_bonds_distances=True)
        # e2 = datetime.datetime.now()-start_t01b
        # m =  "\tTotal time Test01b (Cn with    calc_distances): \t {0:.2f} seconds\n".\
        #     format(float(e2.total_seconds()))
        # m1 = "\t"+len(m)*":"+"\n"
        # print(m1+m+m1) if self.log is None else self.log.info(m1+m+m1)
        # os.rename("test_radiusofgyration/test01/Cn.dat", "test_radiusofgyration/test01/Cn_withdist.dat")
        #
        # # Calculations Cn not calc_distances and Cn with unitvectors ----> Test 01c ==========
        # objcalc_psf = pag.Chain_Statistics(self.trj_small_psf, dt=0.1,
        #                                stride=self.stride, log=self.log)
        # start_t01c = datetime.datetime.now()
        # objcalc = pag.Chain_Statistics(self.trj_small_psf, dt=0.1,
        #                                stride=self.stride, log=self.log)
        # objcalc.calculate(listend2end, diroutput="test_radiusofgyration/test01/",
        #                   acfE2E=True, backbone_list_atoms=backbone_list_atoms,
        #                   isbbatom=is_bb_atoms, calc_Cn_bonds_distances=False,
        #                   single_Cn_unitvector=True)
        # e1 = datetime.datetime.now()-start_t01c
        # m =  "\tTotal time Test01c (Cn without calc_distances and unitvectors): \t {0:.2f} seconds\n".\
        #     format(float(e1.total_seconds()))
        # m1 = "\t"+len(m)*":"+"\n"
        # print(m1+m+m1) if self.log is None else self.log.info(m1+m+m1)
        # os.rename("test_radiusofgyration/test01/Cn.dat", "test_radiusofgyration/test01/Cn_withoutdist_unitvectors.dat")
        #
        # success = self.aux_comparedimensions("test_radiusofgyration/test01/", "RgRee_test01_reference.dat")
        # if success:
        #     m =  "\tComparing Chain dimensions: {}\n".format(success)
        #     m += "\tCalculated values and reference data are equals\n"
        # else:
        #     m = "\tComparing Chain dimensions: {}\n".format(success)
        #     m += "\tWARNING: Calculated vaules and reference data are different\n"
        # print(m) if self.log is None else self.log.info(m)
        #
        # success2 = self.aux_compareinternaldistances("test_radiusofgyration/test01/", "R_kremer_reference.dat")
        # if success2:
        #     m =  "\tComparing Internal chain dimensions: {}\n".format(success2)
        #     m += "\tCalculated values and reference data are equals\n"
        # else:
        #     m = "\tComparing Internal chain dimensions: {}\n".format(success2)
        #     m += "\tWARNING: Calculated vaules and reference data are different\n"
        # print(m) if self.log is None else self.log.info(m)
        #
        # end_t01 = datetime.datetime.now()
        # e = end_t01 - start_t01
        # m =  "\tTotal time Test01: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        # m1 = "\t"+len(m)*"$"+"\n"
        # m2 = "\t============== END   TEST_01 ================================"
        # print(m1+m+m1+m2) if self.log is None else self.log.info(m1+m+m1+m2)

    # ########################################################################
    def test_02_chainstatistics_big(self):

        """
        Test 02: Big chains
        """

        # Create the logger
        start_t02 = datetime.datetime.now()
        m = "\n\t============== START TEST_02 ================================\n"+ \
            "                  Big Trajectory Rg2 and Ree2 \n"
        print(m) if self.log is None else self.log.info(m)

        xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
        xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
        xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
        # xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/trajout.xtc"
        # xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/trajout.xtc"
        # xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/trajout.xtc"

        filename_tpr = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
        filename_psf = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/00-GENERATE/namd_out.psf"
        stride = 1

        # Setup Trajectory
        trj = topology.ExtTrajectory([xtc1, xtc2, xtc3], logger=self.log)
        trj_psf = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

        # Setup topology
        top = topology.Topology(logger=self.log)
        top.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
        trj.set_top_universe(top)

        # Create a list with atoms to calculate the end-to-end distances
        end2end_atoms = []
        nmols_array, neigbors = top.get_array_mols_neigh()
        natoms = trj.topology.natoms
        nchains = len(nmols_array)
        listend2end = []
        natch = 500
        natbb = 500
        for ich in range(nchains):
            ihead = 0 + ich*natch
            itail = (natch-1) + (ich*natch)
            listend2end.append([ich, ihead, itail])

        # Backbone atoms. In this trajectory the branches are located at the end of the chain
        backbone_list_atoms = []
        is_bb_atoms = natoms * [False]
        for ich in range(nchains):
            ihead = 0 + ich * natch
            l = []
            for ibb in range(ihead, ihead+natbb):
                l.append(ibb)
                is_bb_atoms[ibb] = True
            backbone_list_atoms.append(l)

        # Calculations Cn not calc_distances
        start_t01a = datetime.datetime.now()
        objcalc = pag.Chain_Statistics(trj, dt=20000,
                                       stride=stride, log=self.log)
        objcalc.calculate(diroutput="test_radiusofgyration/test02/",
                          listendtoend=listend2end, acfE2E=True, backbone_list_atoms=backbone_list_atoms,
                          isbbatom=is_bb_atoms, calc_Cn_bonds_distances=False)
        e1 = datetime.datetime.now()-start_t01a
        m =  "\tTotal time Test02a (Cn without calc_distances): \t {0:.2f} seconds\n".\
            format(float(e1.total_seconds()))
        m1 = "\t"+len(m)*":"+"\n"
        print(m1+m+m1) if self.log is None else self.log.info(m1+m+m1)
        os.rename("test_radiusofgyration/test02/Cn.dat", "test_radiusofgyration/test02/Cn_withoutdist.dat")

        # Calculations Cn calc_distances
        start_t01b = datetime.datetime.now()
        objcalc = pag.Chain_Statistics(trj, dt=20000,
                                       stride=stride, log=self.log)
        objcalc.calculate(diroutput="test_radiusofgyration/test02/",
                          listendtoend=listend2end, acfE2E=True, backbone_list_atoms=backbone_list_atoms,
                          isbbatom=is_bb_atoms, calc_Cn_bonds_distances=True)
        e2 = datetime.datetime.now()-start_t01b
        m =  "\n\tTotal time Test02b (Cn with    calc_distances): \t {0:.2f} seconds\n".\
            format(float(e2.total_seconds()))
        m1 = "\t"+len(m)*":"+"\n"
        print(m1+m+m1) if self.log is None else self.log.info(m1+m+m1)
        os.rename("test_radiusofgyration/test02/Cn.dat", "test_radiusofgyration/test02/Cn_withdist.dat")

        # Calculations Cn not calc_distances and Cn with unitvectors
        objcalc_psf = pag.Chain_Statistics(trj_psf, dt=20000,
                                       stride=stride, log=self.log)
        start_t01c = datetime.datetime.now()
        objcalc = pag.Chain_Statistics(trj_psf, dt=20000,
                                       stride=stride, log=self.log)
        objcalc.calculate(listend2end, diroutput="test_radiusofgyration/test02/",
                          acfE2E=True, backbone_list_atoms=backbone_list_atoms,
                          isbbatom=is_bb_atoms, calc_Cn_bonds_distances=False,
                          single_Cn_unitvector=True)
        e1 = datetime.datetime.now()-start_t01c
        m =  "\n\tTotal time Test02c (Cn without calc_distances and unitvectors): \t {0:.2f} seconds\n".\
            format(float(e1.total_seconds()))
        m1 = "\t"+len(m)*":"+"\n"
        print(m1+m+m1) if self.log is None else self.log.info(m1+m+m1)
        os.rename("test_radiusofgyration/test02/Cn.dat", "test_radiusofgyration/test02/Cn_withoutdist_unitvectors.dat")

        success = self.aux_comparedimensions("test_radiusofgyration/test02/",
                                             "RgRee_test02_reference.dat", delta=0.1)
        if success:
            m =  "\tComparing Chain dimensions: {}\n".format(success)
            m += "\tCalculated values and reference data are equals\n"
        else:
            m = "\tComparing Chain dimensions: {}\n".format(success)
            m += "\tWARNING: Calculated values and reference data are different\n"
        print(m) if self.log is None else self.log.info(m)

        end_t02 = datetime.datetime.now()
        e = end_t02 - start_t02
        m =  "\tTotal time Test02: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m1 = "\t"+len(m)*"$"+"\n"
        m2 = "\t============== END   TEST_02 ================================"
        print(m1+m+m1+m2) if self.log is None else self.log.info(m1+m+m1+m2)


    # # #########################################################################
    def test_03_chainstatistics_Ree_openmp_single(self):


        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_03 ================================\n"+ \
            "                  OpenMP vs serial C funcitons Ree \n"
        print(m) if self.log is None else self.log.info(m)

        xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
        xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
        xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
        filename_tpr = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
        filename_psf = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/00-GENERATE/namd_out.psf"
        stride = 1

        # Setup Trajectory
        trj = topology.ExtTrajectory([xtc1, xtc2, xtc3], logger=self.log)
        trj_psf = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

        # Setup topology
        top = topology.Topology(logger=self.log)
        top.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
        trj.set_top_universe(top)

        natch = 500
        # Create a list with atoms to calculate the end-to-end distances
        end2end_atoms = []
        nmols_array, neigbors = top.get_array_mols_neigh()
        natoms = trj.topology.natoms
        nchains = len(nmols_array)
        listend2end = []
        natch = 500
        natbb = 500
        for ich in range(nchains):
            ihead = 0 + ich*natch
            itail = (natch-1) + (ich*natch)
            listend2end.append([ich, ihead, itail])


        # Calculations Cn calc_distances
        start_t01b = datetime.datetime.now()
        objcalc = pag.Chain_Statistics(trj, dt=20000,
                                       stride=stride, log=self.log)

        for _ in range(10):
            objcalc.calculate(listend2end, isree=True, isrg=False,
                              molecularweight=True, diroutput="test_radiusofgyration/test03/",
                              unwrap_pbc=True, acfE2E=False, backbone_list_atoms=None,
                              isbbatom=None, calc_Cn_bonds_distances=False,
                              single_Cn_unitvector=False)

        end_t02 = datetime.datetime.now()
        e = end_t02 - start_t01
        m =  "\n\tTotal time Test03: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_03 ================================"
        print(m) if self.log is None else self.log.info(m)

    # # #########################################################################
    def test_04_chainstatistics_Rg_openmp_single(self):

        # print("test_04_chainstatistics_big revisar!!!!!!!!!!!!!!!!!!")
        # return

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_04 ================================\n"+ \
            "                  OpenMP vs serial C functions \n"
        print(m) if self.log is None else self.log.info(m)

        xtc1 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
        xtc2 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
        xtc3 = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"
        filename_tpr = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
        filename_psf = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/00-GENERATE/namd_out.psf"
        stride = 1

        # Setup Trajectory
        trj = topology.ExtTrajectory([xtc1, xtc2, xtc3], logger=self.log)
        trj_psf = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

        # Setup topology
        top = topology.Topology(logger=self.log)
        top.get_bonds_topologyMDAnalysis(filename_tpr, assign_bo=False)
        trj.set_top_universe(top)

        # Calculations Cn calc_distances
        start_t01b = datetime.datetime.now()
        objcalc = pag.Chain_Statistics(trj, dt=20000,
                                       stride=stride, log=self.log)
        for _ in range(10):
            objcalc.calculate(None, molecularweight=True, diroutput="test_radiusofgyration/test04/",
                              unwrap_pbc=True, acfE2E=False, backbone_list_atoms=None,
                              isbbatom=None, calc_Cn_bonds_distances=False,
                              single_Cn_unitvector=False)

    # # #########################################################################
    def test_05_chainstatistics_onesinglechain_nopbc(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_05 ================================\n" + \
            "         Test one single chain without PBC conditions \n"
        print(m) if self.log is None else self.log.info(m)

        gro1 = "../data/P4HB_1chain_pdb_cutoff_verlet/1_4HB_100_box.gro"
        filename_tpr = "../data/P4HB_1chain_pdb_cutoff_verlet/topol.tpr"

        trj_tpr = topology.ExtTrajectory([gro1], topfile=filename_tpr, logger=self.log)
        objcalc = pag.Chain_Statistics(trj_tpr, dt=10, stride=1, log=self.log)

        objcalc.calculate([], diroutput="./", isree=True, isrg=True, iscn=True, acfE2E=False,
                          distributions=False, molecularweight=True, calc_Cn_bonds_distances=True,
                          single_Cn_unitvector=False, begin=0, unwrap_pbc=True,
                          backbone_list_atoms=[], isbondorientation=[],
                          isbbatom=[])

    # ##################################################################################################################
    @classmethod
    def tearDownClass(cls):

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\n\tFinishing at: \t {}".format(now))
        ellapse =  datetime.datetime.now() - cls.start
        cls.log.info("\tTotal time: \t {0:.2f} seconds".format(float(ellapse.total_seconds())))
        m ="\t***************** END RadiusofGyration TEST *****************\n"
        print(m) if cls.log is None else cls.log.info(m)
        utils.close_logger(cls.log)

if __name__ == '__main__':
    unittest.main ()

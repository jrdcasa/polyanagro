import unittest
import datetime
import io
import os
import numpy as np
import polyanagro as pag
import utils
import topology
# from polyanagro.Trajectory import Trajectory
# from polyanagro.Chain_Statistics import Chain_Statistics
# from polyanagro.Topology import Topology

class RDFTests (unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        # Initialize log
        cls.log = utils.init_logger("Output", fileoutput="RDF_test/test_RDF.log",
                                   append=False, inscreen=False)
        cls.log.info("\n\tJob Information\n\t---------------")
        m = "\n\t\t***************** START RDF TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        cls.start = datetime.datetime.now()
        cls.log.info("\t\tStarting: \t {}\n".format(cls.start.strftime("%d-%m-%Y %H:%M:%S")))

    # #########################################################################
    def aux_smalltrj(self, ispsf):

        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
        filename_psf = "../data/0003Ch-C020-002br04/namd_out.psf"

        self.stride = 1

        # Setup Trajectory
        if not ispsf:
            self.trj_small = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_tpr, logger=self.log)
        else:
            self.trj_small = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

    # #########################################################################
    def test_01_RDF_small_bb_tpr(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_01 ================================\n"+ \
            "         Small Trajectory RDF with defaults and tpr topology \n" + \
            "         RDF cannot be performed because bb atoms cannot be setup"
        print(m) if self.log is None else self.log.info(m)

        # Load small trajectory
        self.aux_smalltrj(ispsf=False)

        # Create an object with tpr as topology
        # Create the RDF object. All atoms only backbone in set A and B
        rdfBBAll = pag.RDF(self.trj_small, logger=self.log)
        self.assertFalse(rdfBBAll._isok)

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test01: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_02_RDF_small_bb_psf(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_02 ================================\n"+ \
            "         Small Trajectory RDF with defaults and psf topology \n"
        print(m) if self.log is None else self.log.info(m)

        # Load small trajectory
        self.aux_smalltrj(ispsf=True)

        # Create an object with tpr as topology
        # Create the RDF object. All atoms only backbone in set A and B
        rdfBBAll = pag.RDF(self.trj_small, logger=self.log)
        self.assertTrue(rdfBBAll._isok)
        natomsA = len(rdfBBAll._setA)
        natomsB = len(rdfBBAll._setB)
        self.assertEqual(natomsA, 60)
        self.assertEqual(natomsB, 60)

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test02: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_02 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_03_RDF_small_bbbr_psf(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_03 ================================\n"+ \
            "         Small Trajectory RDF. All atoms and psf topology \n"
        print(m) if self.log is None else self.log.info(m)

        # Load small trajectory
        self.aux_smalltrj(ispsf=True)

        # Create an object with tpr as topology
        # Create the RDF object. All atoms only backbone in set A and B
        rdfBBAll = pag.RDF(self.trj_small, onlybackbone=False, logger=self.log)
        self.assertTrue(rdfBBAll._isok)
        natomsA = len(rdfBBAll._setA)
        natomsB = len(rdfBBAll._setB)
        self.assertEqual(natomsA, 72)
        self.assertEqual(natomsB, 72)

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test02: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_03 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_04_RDF_small_bb_setB_psf(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_04 ================================\n"+ \
            "         Small Trajectory RDF. All atoms and psf topology \n"
        print(m) if self.log is None else self.log.info(m)

        # Load small trajectory
        self.aux_smalltrj(ispsf=True)

        # Create an object with tpr as topology
        # Create the RDF object. Two groups
        rdfBBAll = pag.RDF(self.trj_small, setA_label="Full", setB_label="index 0-10", logger=self.log)
        natomsA = len(rdfBBAll._setA)
        natomsB = len(rdfBBAll._setB)
        self.assertEqual(natomsA, 60)
        self.assertEqual(natomsB, 11)

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test02: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_04 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_05_RDF_small_bb_setAB_psf(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_05 ================================\n"+ \
            "         Small Trajectory RDF. All atoms and psf topology \n"
        print(m) if self.log is None else self.log.info(m)

        # Load small trajectory
        self.aux_smalltrj(ispsf=True)

        # Create an object with tpr as topology
        # Create the RDF object. Two groups
        rdfBBAll = pag.RDF(self.trj_small, setA_label="index 0 1 3 4", setB_label="index 0-10", logger=self.log)
        natomsA = len(rdfBBAll._setA)
        natomsB = len(rdfBBAll._setB)
        self.assertEqual(natomsA, 4)
        self.assertEqual(natomsB, 11)

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test02: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_05 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_06_RDF_small_all_psf_rdf(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_06 ================================\n"+ \
            "         Small Trajectory RDF. All atoms and psf topology \n" +\
            "                    calculate the RDF function"
        print(m) if self.log is None else self.log.info(m)

        # Load small trajectory
        self.aux_smalltrj(ispsf=True)

        # Create an object with tpr as topology
        # Create the RDF object. All atoms all atoms in set A and B
        rdfBBAll = pag.RDF(self.trj_small, onlybackbone=False, logger=self.log)
        rdfBBAll.rdf_calc_cython()

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test02: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_06 ================================"
        print(m) if self.log is None else self.log.info(m)

    # ##################################################################################################################
    @classmethod
    def tearDownClass(cls):

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\n\t\tFinishing at: \t {}".format(now))
        ellapse =  datetime.datetime.now() - cls.start
        cls.log.info("\t\tTotal time: \t {0:.2f} seconds".format(float(ellapse.total_seconds())))
        m ="\t\t***************** END RDF TEST *****************\n"
        print(m) if cls.log is None else cls.log.info(m)
        utils.close_logger(cls.log)

if __name__ == '__main__':
    unittest.main ()

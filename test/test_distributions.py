import unittest
import datetime
import io
import os
import numpy as np
import polyanagro as pag
import topology as top
import utils
import topology


class DistributionTests (unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        # Initialize log
        cls.log = utils.init_logger("Output", fileoutput="Distribution_test/test_distribution.log",
                                   append=False, inscreen=False)
        cls.log.info("\n\tJob Information\n\t---------------")
        m = "\n\t\t***************** START Distribution TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        cls.start = datetime.datetime.now()
        cls.log.info("\t\tStarting: \t {}\n".format(cls.start.strftime("%d-%m-%Y %H:%M:%S")))

    # # #########################################################################
    # def aux_smalltrj(self, ispsf):
    #
    #     xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
    #     xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
    #     xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
    #     filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
    #     filename_psf = "../data/0003Ch-C020-002br04/namd_out.psf"
    #
    #     self.stride = 1
    #
    #     # Setup Trajectory
    #     if not ispsf:
    #         self.trj_small = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_tpr, logger=self.log)
    #     else:
    #         self.trj_small = topology.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=self.log)

    # #########################################################################
    def test_01_Dist_small(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_01 ================================\n"+ \
            "              Create an empty BondedDistributions object  "
        print(m) if self.log is None else self.log.info(m)

        b = pag.BondedDistributions()
        self.assertIsInstance(b, pag.BondedDistributions)

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test01: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        m += "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)


    # #########################################################################
    def test_02_Dist_Bonded(self):

        start_t01 = datetime.datetime.now()
        m = "\n\t============== START TEST_01 ================================\n"+ \
            "              Calculate a bonded distribution       "
        print(m) if self.log is None else self.log.info(m)

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        m = "\t\t Reading topology ({})".format(now)
        print(m) if self.log is None else self.log.info(m)

        filepsf = "./test_data_dist01/namd_out.psf"
        t = top.Topology(logger=self.log)
        t.get_bonds_topologyPSF(filepsf)

        # Arrays for coordinates
        X1 = np.zeros([t.natoms], float)
        Y1 = np.zeros([t.natoms], float)
        Z1 = np.zeros([t.natoms], float)

        b = pag.BondedDistributions()


        #b.bondDist()

        end_t01 = datetime.datetime.now()
        e = end_t01 -start_t01
        m =  "\n\tTotal time Test01: \t {0:.2f} seconds\n".format(float(e.total_seconds()))
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tFinishing: \t {}\n".format(now))
        m += "\t============== END   TEST_01 ================================"
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

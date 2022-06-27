import unittest
import polyanagro as pag
from polyanagro.StemGroup import StemGroup
import topology as top
import utils
import datetime


class StemGroupTests(unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        cls.log = utils.init_logger("Output", fileoutput="StemGroup_test/test_stemgroup.log", append=False, inscreen=False)
        cls.log.info("\n\tJob Information\n\t---------------")
        m = "\n\t\t***************** START Stem Group TEST *****************"
        print(m) if cls.log is None else cls.log.info(m)
        cls.start = datetime.datetime.now()
        now = cls.start.strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\t\tStarting: \t {}\n".format(now))

    # #########################################################################
    def test_01_stemgroup_linear(self):

        m = "\n\t============== START TEST_01 ================================\n"+ \
            "              Test StemGroup with a linear polymer\n"
        print(m) if self.log is None else self.log.info(m)

        dcd1 = "../data/0003Ch-C191-000br00_crystal/traj.dcd"
        filename_psf = "../data/0003Ch-C191-000br00_crystal/namd.psf"
        trj_small_psf = pag.ExtTrajectory([dcd1], topfile=filename_psf, logger=self.log)
        stemObj = StemGroup(trj_small_psf)
        stemObj.look_for_stems()

        m = "\t============== END   TEST_01 ================================"
        print(m) if self.log is None else self.log.info(m)

    # #########################################################################
    def test_02_stemgroup_branch(self):

        m = "\n\t============== START TEST_02 ================================\n"+ \
            "              Test StemGroup with a branch polymer\n"
        print(m) if self.log is None else self.log.info(m)

        dcd1 = "../data/0003Ch-C191-001br04_crystal/traj.dcd"
        filename_psf = "../data/0003Ch-C191-001br04_crystal/namd.psf"
        tempfile = "../data/0003Ch-C191-001br04_crystal/template.dat"
        trj_small_psf = top.ExtTrajectory([dcd1], topfile=filename_psf, templatefile=tempfile, logger=self.log)
        stemObj = StemGroup(trj_small_psf)
        stemObj.look_for_stems()

        m = "\t============== END   TEST_02 ================================"
        print(m) if self.log is None else self.log.info(m)

    # ##################################################################################################################
    @classmethod
    def tearDownClass(cls):

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\n\t\tFinishing at: \t {}".format(now))
        ellapse =  datetime.datetime.now() - cls.start
        cls.log.info("\t\tTotal time: \t {0:.2f} seconds".format(float(ellapse.total_seconds())))
        m ="\t\t***************** END StemGroup TEST *****************\n"
        print(m) if cls.log is None else cls.log.info(m)
        utils.close_logger(cls.log)

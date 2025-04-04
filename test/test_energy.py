import unittest
import datetime
import utils
from polyanagro.Energy import energy_analysis


class EnergyMDTests(unittest.TestCase):

    # #########################################################################
    @classmethod
    def setUpClass(cls):

        # Initialize log
        cls.log = utils.init_logger("Output", fileoutput="test_radiusofgyration/test_chainstatistics.log",
                                   append=False, inscreen=False)
        cls.log.info("\n\tJob Information\n\t---------------")
        m = "\n\t\t***************** START Energy TEST *****************\n"
        m += "\t\t\tEstimated time to run the test is ??? minutes."
        print(m) if cls.log is None else cls.log.info(m)
        cls.start = datetime.datetime.now()
        cls.log.info("\t\tStarting: \t {}\n".format(cls.start.strftime("%d-%m-%Y %H:%M:%S")))


    # ##################################################################################################################
    def test_01_readgromacsenergy(self):

        """
        Read gromacs energy (edr) and make simple plots.

        """

        m = "\t\tTest_01: Check read EDR from GROMACS\n"
        m+= "\t\tThe timestep in the edr is 1ps."
        print(m) if self.log is None else self.log.info(m)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # Create Energy instance
        filename_edr = "../data/energy/GROMACS_edr/set_data_01/ener.edr"

        # All data in edr file is represented in a png file. The files
        # are stored in ./test01_a_figures
        e = energy_analysis(filename_edr, logger=self.log)

        #e = energy_analysis
        e.read_energy()
        e.plot_energy_time(path_to_save="./test_energy/01-test")

    # ##################################################################################################################
    def test_02_gromacsenergybygroup(self):

        """
        Read gromacs energy (edr) and make plots grouped

        """

        m = "\t\tTest_02 Check read EDR from GROMACS and group plot\n"
        m+= "\t\tThe timestep in the edr is 1ps."
        print(m) if self.log is None else self.log.info(m)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # Create Energy instance
        filename_edr = "../data/energy/GROMACS_edr/set_data_01/ener.edr"

        # All data in edr file is represented in a png file. The files
        # are stored in ./test01_a_figures
        e = energy_analysis(filename_edr, logger=self.log)
        e.read_energy()
        e.plot_energy_group(skip_data = 10, path_to_save="./test_energy/02-test")

    # ##################################################################################################################
    def test_03_readgromacsenergy(self):

        """
        Read gromacs energy (edr) and make simple plots.

        """

        m = "\t\tTest_03: Check read EDR from GROMACS (BIG)\n"
        m+= "\t\tThe timestep in the edr is 20ps."
        print(m) if self.log is None else self.log.info(m)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # Create Energy instance
        filename_edr = "../data/energy/GROMACS_edr/set_data_02/ener.part0001.edr"

        # All data in edr file is represented in a png file. The files
        # are stored in ./test01_a_figures
        e = energy_analysis(filename_edr, logger=self.log)

        #e = energy_analysis
        e.read_energy()
        e.plot_energy_time(path_to_save="./test_energy/03-test")

    # ##################################################################################################################
    def test_04_gromacsenergybygroup(self):

        """
        Read gromacs energy (edr) and make plots grouped

        """

        m = "\t\tTest_02 Check read EDR from GROMACS and group plot (BIG)\n"
        m+= "\t\tThe timestep in the edr is 20ps."
        print(m) if self.log is None else self.log.info(m)
        datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")

        # Create Energy instance
        filename_edr = "../data/energy/GROMACS_edr/set_data_02/ener.part0001.edr"

        # All data in edr file is represented in a png file. The files
        # are stored in ./test01_a_figures
        e = energy_analysis(filename_edr, logger=self.log)
        e.read_energy()
        e.plot_energy_group(skip_data = 10, path_to_save="./test_energy/04-test")

    # ##################################################################################################################
    @classmethod
    def tearDownClass(cls):

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        cls.log.info("\n\tFinishing at: \t {}".format(now))
        ellapse =  datetime.datetime.now() - cls.start
        cls.log.info("\tTotal time: \t {0:.2f} seconds".format(float(ellapse.total_seconds())))
        m ="\t***************** END Energy TEST *****************\n"
        print(m) if cls.log is None else cls.log.info(m)
        utils.close_logger(cls.log)

if __name__ == '__main__':
    unittest.main ()
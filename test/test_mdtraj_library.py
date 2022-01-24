import unittest
import polyanagro as pag
import numpy as np
from copy import copy
import warnings
import datetime

class MDTRAJTests(unittest.TestCase):

    # ##################################################################################################################
    @classmethod
    def setUpClass(self):

        # Initialize log
        self.log = pag.init_logger("Output", fileoutput="trj_test/test_mdtraj_library.log",
                                   append=False, inscreen=False)
        self.log.info("\n\tJob Information\n\t---------------")
        warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
        m = "\n\t\t***************** START MDTRAJ TEST *****************"
        print(m) if self.log is None else self.log.info(m)
        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\t\tStarting: \t {}\n".format(now))

    # # #################################################################################################################
    def test_01_UApolTrappe(self):

        m = "\t============== TEST_01 TrappeUA Polymer ================================\n"+ \
            "\t   Testing to read trajectories of a polymer typed with Trappe-UA\n"+ \
            "\t   Basic test for topology and trajectory read\n"
        print(m) if self.log is None else self.log.info(m)

        # List of trajectory files
        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"

        # Topology files
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
        filename_psf = "../data/0003Ch-C020-002br04/namd.psf"

        # Create a trajectory with all trajectory files
        trj = pag.ExtTrajectory([xtc1, xtc2, xtc3], logger=self.log)

        # Create topology with tpr file
        top1 = pag.Topology(logger=self.log)
        top1.get_bonds_topologyMDAnalysis(filename_tpr)
        top1.draw_graph_forest_networkx("trj_test/test_01_")

        # Create topology with psf file and MDAnalysis library
        top2 = pag.Topology(logger=self.log)
        top2.get_bonds_topologyMDAnalysis(filename_psf)

        # Create topology with psf file
        top3 = pag.Topology(logger=self.log)
        top3.get_bonds_topologyPSF(filename_psf)

        m = "\ttop1 and top2 are equal: {}".format(top1 == top2)
        print(m) if self.log is None else self.log.info(m)

        self.assertNotEqual(top1, top2)
        self.assertNotEqual(top1, top3)
        self.assertEqual(top2, top3)

        self.assertEqual(trj.trjtype, ['xtc','xtc','xtc'])

        top1.assign_bond_orders()
        # There is not order bonds greater than 1.0
        l1 = top1._orderbonds[np.where(top1._orderbonds >1.0)]
        self.assertEqual(list(l1),[])
        # There is any order bonds equals to 1.0
        l2 = top1._orderbonds[np.where(top1._orderbonds ==1.0)]
        self.assertNotEqual(list(l2),[])

        m = "\t==============        END   TEST_01       ================================"
        print(m) if self.log is None else self.log.info(m)

    # #################################################################################################################
    def test_02_topologycopy(self):

        m = "\n\t============== TEST_02 Copy of topology and trajectory ================================\n"+ \
            "\t   Testing the copy of a topology and trajectory\n"
        print(m) if self.log is None else self.log.info(m)

        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"

        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"

        # Trajectory copy
        trj = pag.ExtTrajectory([xtc1, xtc2, xtc3], logger=self.log)
        trj_copy = copy(trj)

        self.assertEqual(trj.trjtype, trj_copy.trjtype)
        self.assertEqual(trj.trjpath, trj_copy.trjpath)
        self.assertEqual(trj.topology, trj_copy.topology)

        # Topology copy
        top1 = pag.Topology(logger=self.log)
        top1.get_bonds_topologyMDAnalysis(filename_tpr)
        top1.draw_graph_forest_networkx("trj_test/test_02_")

        top1_copy = copy(top1)

        self.assertEqual(top1, top1_copy)
        top1_copy._nmols = []
        self.assertNotEqual(top1, top1_copy)

        m = "\t==============        END   TEST_02       ================================"
        print(m) if self.log is None else self.log.info(m)

    #################################################################################################################
    def test_03_OPLSEthanolWater(self):

        m = "\n\t============== TEST_03 Ethanol/Water OPLS ================================\n"+ \
            "\t Testing to read trajectories of a mix of small molecules using OPLS\n"
        print(m) if self.log is None else self.log.info(m)

        xtc1 = "../data/ethanol-water/RUN-0001/traj_comp.xtc"
        filename_tpr = "../data/ethanol-water/RUN-0001/topol.tpr"

        trj = pag.ExtTrajectory([xtc1], logger=self.log)
        top1 = pag.Topology(logger=self.log)
        top1.get_bonds_topologyMDAnalysis(filename_tpr)

        self.assertEqual(trj.trjtype, ['xtc'])

        top1.assign_bond_orders()
        # There is not order bonds greater than 1.0
        l1 = top1._orderbonds[np.where(top1._orderbonds >1.0)]
        self.assertEqual(list(l1),[])
        # There is any order bonds equals to 1.0
        l2 = top1._orderbonds[np.where(top1._orderbonds ==1.0)]
        self.assertNotEqual(list(l2),[])

        m = "\t==============        END   TEST_03       ================================"
        print(m) if self.log is None else self.log.info(m)

    ##################################################################################################################
    def test_04_Aromatic(self):

        m = "\n\t============== TEST_04 Aromatic molecule ================================\n"+ \
            "\t   Testing to read trajectories of containing aromatic molecules (benzene)\n"
        print(m) if self.log is None else self.log.info(m)

        xtc1 = "../data/benzene/RUN-0001/traj_comp.xtc"
        filename_tpr = "../data/benzene/RUN-0001/topol.tpr"

        trj = pag.ExtTrajectory([xtc1], logger=self.log)
        top1 = pag.Topology(logger=self.log)
        top1.get_bonds_topologyMDAnalysis(filename_tpr)
        top1.draw_graph_forest_networkx("trj_test/test_04")

        top1.assign_bond_orders()
        # There is not order bonds greater than 1.0
        l1 = top1._orderbonds[np.where(top1._orderbonds  == 1.5)]
        self.assertNotEqual(list(l1),[])
        # There is any order bonds equals to 1.0
        l2 = top1._orderbonds[np.where(top1._orderbonds ==1.0)]
        self.assertNotEqual(list(l2),[])

        m = "\t==============        END   TEST_04       ================================"
        print(m) if self.log is None else self.log.info(m)

    #################################################################################################################
    def test_05_trj(self):

        m = "\n\t============== TEST_05 Check Read multiple trajectories ================================\n"+ \
            "\t   Check Read Multiple Trajectories \n"
        print(m) if self.log is None else self.log.info(m)

        # Input files
        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"

        # Setup Trajectory
        trj = pag.ExtTrajectory([xtc1, xtc2, xtc3], logger=self.log)

        # Setup topology
        top = pag.Topology(logger=self.log)
        top.get_bonds_topologyMDAnalysis(filename_tpr)
        trj.set_top_universe(top)

        # Write trajectory xtc (by extension)
        trj.write_trajectory("trj_test/traj_test05_jump_w.xtc", nojump=False)
        trj.write_trajectory("trj_test/traj_test05_nojump_w.xtc", nojump=True)

        # Write trajectory dcd (by format)
        trj.write_trajectory("trj_test/traj_test05_unwrap_w.dcd", format="dcd")

        m = "\t==============        END   TEST_05       ================================"
        print(m) if self.log is None else self.log.info(m)

    # #################################################################################################################
    def test_06_trj_with_topology(self):

        m = "\n\t============== TEST_06 Check Read multiple trajectories and topology ================================\n"+ \
            "     Check Read Multiple Trajectories with Topology \n"
        print(m) if self.log is None else self.log.info(m)

        # Input files
        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"

        # Setup Trajectory and Topology
        trj = pag.ExtTrajectory([xtc1, xtc2, xtc3], filename_tpr, logger=self.log)
        # Write trajectory xtc (by extension)
        trj.write_trajectory("trj_test/trajtopo_test06_nojump_w.xtc", nojump=True)

        nframes = trj.get_numframes()
        self.assertEqual(nframes, 301)

        m = "\t==============        END   TEST_06       ================================"
        print(m) if self.log is None else self.log.info(m)

    # # ##################################################################################################################
    def test_07_joinbigxtc_trjs(self):

        r"""Test read XTC trajectories with MDAnalysis

        Time write unwrapped trajectory XTC: 325 s --> 5min 25s (unwrap MDAnalysis function)
        Time write unwrapped trajectory XTC:  27 s --> 0min 27s (unwrap Cython function)
        Time write unwrapped trajectory XTC:  87 s --> 1min 20s (unwrap nojump Cython function)
        Time write   wrapped trajectory XTC:  21 s
        Time write   wrapped trajectory DCD:  19 s

        """
        m = "\n\t============== TEST_07 Check Read multiple big trajectories ================================\n"+ \
            "\t     Check Read Multiple big Trajectories \n"
        print(m) if self.log is None else self.log.info(m)

        tpr_name = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/topol.tpr"
        xtc1_name = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/01-RESTART-0000-1000ns/traj_comp.xtc"
        xtc2_name = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/02-RESTART-1000-2000ns/traj_comp.part0002.xtc"
        xtc3_name = "../data/0100Ch-C500-000br02/T-450K/06-NPT-100ns-Parrinello/03-RESTART-2000-3000ns/traj_comp.part0003.xtc"

        xtc_list = [xtc1_name, xtc2_name, xtc3_name]

        trj = pag.ExtTrajectory(xtc_list, topfile=tpr_name, logger=self.log)

        # Written trajectory is of the same type that input trajectory
        start_time = datetime.datetime.now()
        trj.write_trajectory("trj_test/trjfull_out_big_unwrap.xtc", pbc=True)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\tUnwrapping coordinates --> time: {0:s} seconds".format(str(elapsed_time.total_seconds()))
        self.log.info(m)

        #Written trajectory is of the same type that input trajectory
        start_time = datetime.datetime.now()
        trj.write_trajectory("trj_test/trjfull_out_big_unwrap_nojump.xtc", pbc=True, nojump=True)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\tUnwrapping No Jump coordinates --> time: {0:s} seconds".format(str(elapsed_time.total_seconds()))
        self.log.info(m)

        # Written trajectory is of the same type that input trajectory
        start_time = datetime.datetime.now()
        trj.write_trajectory("trj_test/trjfull_out_big.xtc", pbc=False, nojump=False)
        end_time = datetime.datetime.now()
        elapsed_time = end_time - start_time
        m = "\tNo Wrapping coordinates --> time: {0:s} seconds".format(str(elapsed_time.total_seconds()))
        self.log.info(m)

        m = "\t==============        END   TEST_07       ================================"
        print(m) if self.log is None else self.log.info(m)

    ################################################################################################################
    def test_08_trj_with_topology_addtrj(self):

        m = "\n\t============== TEST_08 Add trajectory ================================\n"+ \
            "     Check add a new trajectory file to a Trajectory with Topology \n"
        print(m) if self.log is None else self.log.info(m)

        # Input files
        xtc1 = "../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
        xtc2 = "../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
        filename_tpr = "../data/0003Ch-C020-002br04/RUN-001/topol.tpr"

        # Setup Trajectory and Topology
        trj = pag.ExtTrajectory([xtc1, xtc2], filename_tpr, logger=self.log)
        # Write trajectory xtc (by extension)
        nframes = trj.get_numframes()
        self.assertEqual(nframes, 201)

        xtc3 = "../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
        trj.add_trjs(xtc3)
        trj.universe.load_new(trj.trjpath, continuous=True)
        nframes = trj.get_numframes()
        self.assertEqual(nframes, 301)

        m = "\t==============        END   TEST_08       ================================"
        print(m) if self.log is None else self.log.info(m)

    # ##################################################################################################################
    def test_99(self):

        now = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S")
        self.log.info("\n\t\tFinishing at: \t {}".format(now))
        m ="\t\t***************** END MDTRAJ TEST *****************\n"
        print(m) if self.log is None else self.log.info(m)

        pag.close_logger(self.log)

    # ##################################################################################################################
    def tearDown(self):

        pass



if __name__ == '__main__':

    unittest.main()

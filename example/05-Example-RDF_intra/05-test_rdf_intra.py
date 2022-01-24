import polyanagro as pag

# Create the logger
# Logger
log = pag.init_logger("Output", fileoutput="test_odf.log",
                          append=False, inscreen=True)
# End Logger

# Input files
xtc1 = "../../data/0003Ch-C020-002br04/RUN-001/traj_comp.xtc"
xtc2 = "../../data/0003Ch-C020-002br04/RUN-002/traj_comp.part0002.xtc"
xtc3 = "../../data/0003Ch-C020-002br04/RUN-003/traj_comp.part0003.xtc"
filename_tpr = "../../data/0003Ch-C020-002br04/RUN-001/topol.tpr"
filename_psf = "../../data/0003Ch-C020-002br04/namd_out.psf"
stride = 1000

# Setup Trajectory
trj = pag.ExtTrajectory([xtc1, xtc2, xtc3], topfile=filename_psf, logger=log)
print(trj.topology._isbackbone)


# Create the RDF object. All atoms only backbone in set A and B
rdfBBAll = pag.RDF(trj, logger=log)
#rdfBBAll.rdf_calc_cython()

# RDF only C2 atoms
#rdfstats.rdf_calc_cython()


print("Job Done!!!!")



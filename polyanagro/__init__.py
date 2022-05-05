"""
`polyanagro` analyzes polymer trajectories from MD or MC calculations
"""
from polyanagro.version import __version__
from polyanagro.Calculations import Calculations
from polyanagro.Chain_Statistics import Chain_Statistics
from polyanagro.RDF import RDF
from polyanagro.BondedDistribution import BondedDistributions

from topology.atomic_data import element_cov_radius, maximal_valences, \
                                   united_atoms_equivalence, atomic_mass, \
                                   element_vdw_vmd_volume_bondi

from ext_libc.c_rg_openmp import calc_rg_openmp_massweigth, calc_rg_openmp
from ext_libc.c_unit_bond_vectors import unit_bond_vectors, bond_bond_correlation
from ext_libc.c_rdf_openmp import setup_rdf_init, rdf_hist, rdf_gr
from ext_libc.c_distC import setup_hist_bondC, setup_hist_angleC, setup_hist_dihC, bondDistC
from ext_libc.c_acf_openmp import calc_acf_ete



#
# from polyanagro.utils import padding_list, total_size, get_name_file_ext
#
# from polyanagro.StemGroup import StemGroup
#
#
#
# from polyanagro.internal_coordinates import distance_array, generate_random_euler_angles,\
#                                             euler_rotation_matrix, center_of_geom_purepython,\
#                                             unwrap_purepython, unwrap, unwrap_nojump, \
#                                             dihedral_angle_purepython, bend_angle_purepython,\
#                                             cos_angle_purepython, distance_array_purepython,\
#                                             distance_array_numpypython, distance_diagonal_array,\
#                                             center_of_geom, center_of_mass_purepython, center_of_mass,\
#                                             calc_distance_array_numba
#
# from polyanagro.internal_coordinates_numba import unwrap_numba
#


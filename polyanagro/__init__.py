"""
`polyanagro` analyzes polymer trajectories from MD or MC calculations
"""
from polyanagro.version import __version__
from polyanagro.Calculations import Calculations
from polyanagro.Chain_Statistics import Chain_Statistics
from polyanagro.RDF import RDF
from polyanagro.BondedDistribution import BondedDistributions
from polyanagro.Neighbors import Neighbors


from topology.atomic_data import element_cov_radius, maximal_valences, \
                                   united_atoms_equivalence, atomic_mass, \
                                   element_vdw_vmd_volume_bondi

from ext_libc.c_rg_openmp import calc_rg_openmp_massweigth, calc_rg_openmp
from ext_libc.c_unit_bond_vectors import unit_bond_vectors, bond_bond_correlation
from ext_libc.c_rdf_openmp import setup_rdf_init, rdf_hist, rdf_gr
from ext_libc.c_distC import setup_hist_bondC, setup_hist_angleC, setup_hist_dihC, \
                             bondDistC, angleDistC, dihDistC
from ext_libc.c_acf_openmp import calc_acf_ete
from ext_libc.c_unit_bond_vectors import setup_odf_intra, odf_intra, avg_write_odf_intra
from ext_libc.c_internal_distances import setup_internal_distances, \
                                          internal_distances_iframe, \
                                          internal_distances_return,\
                                          insternalchaindist_cython




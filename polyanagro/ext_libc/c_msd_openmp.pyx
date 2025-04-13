# import both numpy and the Cython declarations for numpy
import cython
import datetime
from libc.stdlib cimport malloc, free
from time import time
from cython.parallel import prange
import numpy as np
cimport numpy as np
from libc.stdio cimport fopen, fprintf, fclose, FILE
cdef bint USED_OPENMP

# declare the interface to the C code
cdef extern from "calc_msd.c":
    cdef bint USED_OPENMP
    #void c_msd_fft3(double*, double*, int, int)
    void c_msd_all(double*, int, int, int, double, int, const char*)
cdef extern from "calc_msd_fftw3.c":
    double* c_msd_fftw3_fast(double *, int, int, int)

OPENMP_ENABLED = True if USED_OPENMP else False

# ========================================================================================
@cython.boundscheck(False)
@cython.wraparound(False)
def msd_atomistic_cython(int nframes, int natoms, int dim, float dt, float[:, :, :] positions, bytes filename):

    cdef int i, j, irow, icol
    cdef float** disp
    cdef float* sq_disp
    cdef char* c_filename = filename  # Convert bytes to char*

    # Allocate memory ============
    disp = <float**> malloc(natoms * sizeof(float*))
    if disp is NULL:
        raise MemoryError("Failed to allocate memory for disp")

    for i in range(natoms):
        disp[i] = <float*> malloc(dim * sizeof(float))
        if disp[i] is NULL:
            raise MemoryError("Failed to allocate memory for disp[i]")
        for j in range(dim):
            disp[i][j] = 0.0  # Initialize to zero

    sq_disp = <float*> malloc(natoms * sizeof(float))  # Fix allocation
    if sq_disp is NULL:
        raise MemoryError("Failed to allocate memory for sq_disp")
    for i in range(natoms):
        sq_disp[i] = 0.0

    msd = <float*> malloc(nframes * sizeof(float))  # Fix allocation
    if msd is NULL:
        raise MemoryError("Failed to allocate memory for msd")

    count = <int*> malloc(nframes * sizeof(int))  # Fix allocation
    if count is NULL:
        raise MemoryError("Failed to allocate memory for count")

    times = <float*> malloc(nframes * sizeof(float))  # Fix allocation
    if times is NULL:
        raise MemoryError("Failed to allocate memory for times")


    for i in range(nframes):
        msd[i] = 0.0
        count[i] = 0
        times[i] = 0.0

    # MSD Calculation
    start_time = datetime.datetime.now()
    for i in range(0, nframes):

        if i % 100 == 0:
            msg = "\tFrame {0:9d} of {1:9d}".format(i, nframes)
            mid_time = datetime.datetime.now()
            elapsed_time = mid_time - start_time
            msg += "\ttime: {0:s} seconds".format(str(elapsed_time.total_seconds()))
            print(msg)

        for j in range(i, nframes):
            for irow in range(0, natoms):

                for icol in range(0, dim):
                    disp[irow][icol] = positions[j][irow][icol] - positions[i][irow][icol]
                sq_disp[irow] = disp[irow][0]*disp[irow][0] + \
                                disp[irow][1]*disp[irow][1] + \
                                disp[irow][2]*disp[irow][2]
                msd[j-i] += sq_disp[irow]
            count[j-i] += natoms

    for i in range(nframes):
        msd[i] = msd[i] / float(count[i])
        times[i] = i * dt

    # Print
    with open(c_filename, 'w') as f:
        for i in range(nframes):
            line = "{0:.3f}  {1:.6f}\n".format(times[i], msd[i])
            f.writelines(line)

    # Free memory ============
    if disp is not NULL:
        for i in range(natoms):
            if disp[i] is not NULL:
                free(disp[i])  # Free inner arrays
        free(disp)  # Free outer array

    if sq_disp is not NULL:
        free(sq_disp)  # Free 1D array

    if times is not NULL:
        free(times)  # Free 1D array

    if count is not NULL:
        free(count)  # Free 1D array

    if msd is not NULL:  # FIX: Free `msd`
        free(msd)

# ========================================================================================
# @cython.boundscheck(False)
# @cython.wraparound(False)
# @cython.cdivision(True)  # Avoids unnecessary division checks
def msd_atomistic_opt_cython(int nframes, int natoms, int dim, float dt, float[:, :, :] positions, bytes filename):

    cdef int i, j, irow, icol
    cdef float **disp
    cdef float *sq_disp, *msd, *times
    cdef int *count
    cdef char *c_filename = filename  # Convert bytes to char*

    # Allocate memory using numpy (avoids malloc/free issues)
    disp = <float**> malloc(natoms * sizeof(float*))
    for i in range(natoms):
        disp[i] = <float*> malloc(dim * sizeof(float))

    sq_disp = <float*> malloc(natoms * sizeof(float))
    msd = <float*> malloc(nframes * sizeof(float))
    count = <int*> malloc(nframes * sizeof(int))
    times = <float*> malloc(nframes * sizeof(float))

    if not disp or not sq_disp or not msd or not count or not times:
        raise MemoryError("Memory allocation failed")

    # Initialize arrays
    for i in range(nframes):
        msd[i] = 0.0
        count[i] = 0
        times[i] = i * dt

    # Parallelized MSD Calculation
    start_time = time()

    with nogil:
        for i in prange(nframes, schedule="static"):  # OpenMP parallelized
            for j in range(i, nframes):
                for irow in range(natoms):

                    # Compute displacement
                    for icol in range(dim):
                        disp[irow][icol] = positions[j, irow, icol] - positions[i, irow, icol]

                    # Compute squared displacement dynamically (supports any `dim`)
                    sq_disp[irow] = 0.0
                    for icol in range(dim):
                        sq_disp[irow] += disp[irow][icol] * disp[irow][icol]

                    # Accumulate MSD
                    msd[j - i] += sq_disp[irow]

                count[j - i] += natoms

        # Normalize MSD
        for i in range(nframes):
            msd[i] /= float(count[i])

    end_time = time()
    elapsed_time = end_time - start_time
    msg = "\ttime: {0:s} seconds".format(str(elapsed_time))
    print(msg)

    # Write results to file (avoids Python overhead)
    cdef FILE* f = fopen(c_filename, "w")
    if f != NULL:
        for i in range(nframes):
            fprintf(f, "%.3f  %.6f\n", times[i], msd[i])
        fclose(f)

    # Free allocated memory
    for i in range(natoms):
        free(disp[i])
    free(disp)
    free(sq_disp)
    free(msd)
    free(count)
    free(times)


# =============================================================================
def msd_all_c(np.ndarray[np.float64_t, ndim=3] positions,
              int nmol, double timestep, int weedf, str filename):

    """

    Args:
        positions: (nframes, natoms, 3) -> Positions of all atoms in the trayectory
        nmol: Number of molecules
        timestep: Timestep in ps
        weedf: Number of starting points of independient time series
        filename: Name for the MSD filename

    Returns:

    """

    cdef int nframes = positions.shape[0]  # Number of frames
    cdef int natoms = positions.shape[1]  # Number of atoms

    cdef bytes encoded = filename.encode('utf-8')  # keep it alive
    cdef const char *c_str = encoded

    c_msd_all(&positions[0, 0, 0], nframes, nmol, natoms, timestep, weedf, c_str)


# =============================================================================
def msd_fftw3_cython(np.ndarray[np.float64_t, ndim=3, mode='c'] positions,
                     double timestep):

    cdef int natoms = positions.shape[0]  # nSignals
    cdef int ndim = positions.shape[1]    # Dimnesions
    cdef int nframes = positions.shape[2] # Signalsize
    #cdef int nsignals = natoms*ndim
    # cdef int signalsize = nframes
    cdef double* result

    # Check memory use

    result = c_msd_fftw3_fast(&positions[0, 0, 0], natoms, nframes , ndim)

    if result == NULL:
        raise MemoryError("computeMSD returned NULL")

    # Total length of result array
    cdef int size = nframes * ndim
    # Wrap into NumPy array (copy it)
    msd = np.empty(size, dtype=np.float64)
    for i in range(size):
        msd[i] = result[i]

    return msd

# =============================================================================
def msd_com_c(np.ndarray[int, ndim=2, mode="c"] mols,
              np.ndarray[np.float64_t, ndim=1, mode="c"] mass,
              np.ndarray[np.float64_t, ndim=3, mode="c"] positions_atoms,
              np.ndarray[np.float64_t, ndim=3, mode="c"] positions_com):

    cdef int nchains = mols.shape[0]
    cdef int maxatomsch = mols.shape[1]
    cdef int nframes = positions_atoms.shape[0]
    cdef int dimensions = positions_atoms.shape[2]

    cdef int ich, i, iframe, d
    cdef int atom_index
    cdef double m, total_mass
    cdef double tmp_pos[3]

    for iframe in range(nframes):
        for ich in range(nchains):
            total_mass = 0.0
            tmp_pos[0] = tmp_pos[1] = tmp_pos[2] = 0.0

            for i in range(maxatomsch):
                atom_index = mols[ich, i]
                if atom_index < 0:
                    continue  # skip invalid entries

                m = mass[atom_index]
                total_mass += m
                for d in range(dimensions):
                    tmp_pos[d] += m * positions_atoms[iframe, atom_index, d]

            if total_mass > 0.0:
                for d in range(dimensions):
                    positions_com[iframe, ich, d] = tmp_pos[d] / total_mass
            else:
                for d in range(dimensions):
                    positions_com[iframe, ich, d] = 0.0  # fallback if mass sum is zero

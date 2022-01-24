import numpy as np
import numba
from numba import jit
import time

SQRT_2PI = np.sqrt(2*np.pi)

@jit(nopython=True, parallel=True, cache=True)
def gaussians(x, means, widths):

    """
    Return the value of gaussian kernels.

        x - location of evaluation
        means - array of kernel means
        widths - array of kernel widths
    """
    n = means.shape[0]
    result = np.exp(-0.5 * ((x - means) / widths)**2) / widths

    return result/SQRT_2PI/n

@jit(nopython=True, cache=True)
def gaussians_nothread(x, means, widths):

    """
    Return the value of gaussian kernels.

        x - location of evaluation
        means - array of kernel means
        widths - array of kernel widths
    """
    n = means.shape[0]
    result = np.exp(-0.5 * ((x - means) / widths)**2) / widths

    return result/SQRT_2PI/n

def gaussians_numpy(x, means, widths):

    """
    Return the value of gaussian kernels.

        x - location of evaluation
        means - array of kernel means
        widths - array of kernel widths
    """
    n = means.shape[0]
    result = np.exp(-0.5 * ((x - means) / widths)**2) / widths

    return result/SQRT_2PI/n

if __name__ == "__main__":

    """
    Results for 100 iterations with cache=False ==============
        Numpy pure = 1.876 seconds (S = 1.00)
        No Thread  = 1.399 seconds (S = 1.34)
        Thread     = 0.552 seconds (S = 3.40)
    Results for 100 iterations with cache=True ============== 
        Numpy pure = 1.862 seconds (S = 1.00)
        No Thread  = 1.132 seconds (S = 1.64)
        Thread     = 0.225 seconds (S = 8.27)    

    Results for 1000 iterations with cache=False ===============
        Numpy pure = 16.769 seconds (S = 1.00)
        No Thread  = 10.365 seconds (S = 1.62)
        Thread     = 2.407  seconds (S = 6.97)
 
    Results for 1000 iterations with cache=True ============== 
        Numpy pure = 18.491 seconds (S = 1.00)
        No Thread  = 10.060 seconds (S = 1.84)
        Thread     = 2.156  seconds (S = 8.58)   
    
    Results for 5000 iterations with cache=False ==============
        Numpy pure = 91.298 seconds (S = 1.00)
        No Thread  = 49.378 seconds (S = 1.85)
        Thread     = 11.280 seconds (S = 8.09)
    Results for 5000 iterations with cache=True ============== 
        Numpy pure = 87.886 seconds (S = 1.00)
        No Thread  = 48.873 seconds (S = 1.80)
        Thread     = 10.582 seconds (S = 8.31)

    """

    means = np.random.uniform(-1, 1, size=1000000)
    widths = np.random.uniform(0.1, 0.3, size=1000000)

    iter = 1000

    print("Results for {0:d} iterations with cache=False ".format(iter))

    ts3 = time.time()
    for _ in range(iter):
        r3 = gaussians_numpy(0.4, means, widths)
    tf3 = time.time() - ts3
    print("\tNumpy pure = {0:.3f} seconds (S = {1:.2f})".format(tf3,tf3/tf3))

    ts2 = time.time()
    #gaussians_nothread =jit(nopython=True, cache=True, parallel=False)(gaussians.py_func)
    for _ in range(iter):
        r2 = gaussians_nothread(0.4, means, widths)
    tf2 = time.time() - ts2
    print("\tNo Thread  = {0:.3f} seconds (S = {1:.2f})".format(tf2,tf3/tf2))

    ts1 = time.time()
    for _ in range(iter):
        r1 = gaussians(0.4, means, widths)
    tf1 = time.time() - ts1
    print("\tThread     = {0:.3f} seconds (S = {1:.2f})".format(tf1,tf3/tf1))





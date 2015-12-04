import numpy as np
cimport numpy as np
import cython

cimport cmicroscope

@cython.boundscheck(False)
@cython.wraparound(False)
def emission(
        np.ndarray[double, ndim=2, mode="c"] points not None,
        np.ndarray[double, ndim=1, mode="c"] intensity not None,
        data_size = None,
        double ATsq = 1.16737230263e+25, double d = 1.01896194663e+03,
        double epsilon = 1.0, double absorption_cross_section = 2e-18,
        double exposure_time = 100e-3):
    if data_size is None:
        data_size = intensity.size
    return cmicroscope.emission(
        <double (*)[3]>&points[0, 0], &intensity[0], <unsigned int>data_size,
        ATsq, d, epsilon, absorption_cross_section, exposure_time)

@cython.boundscheck(False)
@cython.wraparound(False)
def overlay_psf(
        np.ndarray[double, ndim=1, mode="c"] data not None,
        unsigned int N_pixel, double pixel_length,
        np.ndarray[double, ndim=1, mode="c"] p not None,
        double I,
        np.ndarray[double, ndim=1, mode="c"] c not None,
        double k, double N_A, double cutoff):
    cmicroscope.overlay_psf(&data[0], N_pixel, pixel_length, &p[0], I, &c[0], k, N_A, cutoff)

@cython.boundscheck(False)
@cython.wraparound(False)
def detection(
        np.ndarray[double, ndim=1, mode="c"] input not None,
        np.ndarray[double, ndim=1, mode="c"] output not None,
        data_size = None):
    if data_size is None:
        data_size = input.size
    cmicroscope.detection(&input[0], &output[0], <unsigned int>data_size)

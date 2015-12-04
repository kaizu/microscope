cdef extern from "microscope.hpp" namespace "microscope":
    double emission(double[][3], double[], unsigned int, double, double, double, double, double)
    void overlay_psf(double[], unsigned int, double, double[3], double, double[3], double, double, double)
    void detection(double[], double[], unsigned int)

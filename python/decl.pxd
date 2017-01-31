from libcpp.pair cimport pair

cdef extern from "microscope.hpp" namespace "microscope":
    pair[double, double] tirf_evanescent_field(double, double, double)
    double emission(double[][3], double[], unsigned int, double, double, double, double, double, double)
    void overlay_psf(double[], unsigned int, double, double[3], double, double[3], double, double, double)
    void cmos_detection(double[], double[], unsigned int)
    void emccd_detection(double[], double[], unsigned int)

    void detection(double[], double[], unsigned int) #DEPRECATED: == cmos_detection

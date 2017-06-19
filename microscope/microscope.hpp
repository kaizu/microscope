#ifndef __MICROSCOPE__MICROSCOPE
#define __MICROSCOPE__MICROSCOPE

#include <algorithm>
#include <gsl/gsl_rng.h>

#include "point_spreading_functions.hpp"
#include "emission.hpp"
#include "detection.hpp"

#include <time.h>

namespace microscope
{

double emission(
    double points[][3], double intensity[], unsigned int const data_size,
    double const ATsq = 1.16737230263e+25, double const d = 1.01896194663e+03,
    double const epsilon = 1.0, double const absorption_cross_section = 2e-18,
    double const exposure_time = 100e-3, double const QY = 0.61)
{
    double Itot = 0.0;
    for (unsigned int i(0); i < data_size; ++i)
    {
        double const z = fabs(points[i][2]);
        double const I = tirf_emission(
            z, ATsq, d, epsilon, absorption_cross_section, exposure_time);
        intensity[i] *= I;
        Itot += intensity[i];
    }
    return Itot;
}

void overlay_psf(
    double data[], unsigned int const N_pixel, double const pixel_length,
    double p[3], double const I, double c[3], double const k, double const N_A,
    double const cutoff)
{
    const double offset(N_pixel * pixel_length * -0.5);
    const double x(p[0] - c[0] - offset);
    const double y(p[1] - c[1] - offset);
    const unsigned int imin(
        static_cast<unsigned int>(
            std::max(0.0, floor((x - cutoff) / pixel_length))));
    const unsigned int imax(
        std::min(N_pixel,
            static_cast<unsigned int>(ceil((x + cutoff) / pixel_length))));
    const unsigned int jmin(
        static_cast<unsigned int>(
            std::max(0.0, floor((y - cutoff) / pixel_length))));
    const unsigned int jmax(
        std::min(N_pixel,
            static_cast<unsigned int>(ceil((y + cutoff) / pixel_length))));

    for (unsigned int i(imin); i < imax; ++i)
    {
        const double xmin(pixel_length * i + offset);
        const double xmax(xmin + pixel_length);

        for (unsigned int j(jmin); j < jmax; ++j)
        {
            const double ymin(pixel_length * j + offset);
            const double ymax(ymin + pixel_length);

            // const double value(int_psf_gaussian(
            //     xmin, xmax, ymin, ymax, p, c, k, N_A));
            const double value(int_psf_tbl(
                xmin, xmax, ymin, ymax, p, c, k, N_A));
            data[i * N_pixel + j] += I * value;
        }
    }

    // for (unsigned int i(0); i < N_pixel; ++i)
    // {
    //     const double xmin(pixel_length * (i - N_pixel * 0.5));
    //     const double xmax(xmin + pixel_length);
    //     for (unsigned int j(0); j < N_pixel; ++j)
    //     {
    //         const double ymin(pixel_length * (j - N_pixel * 0.5));
    //         const double ymax(ymin + pixel_length);

    //         const double value(int_psf_with_cutoff(
    //             p, c, k, N_A, xmin, xmax, ymin, ymax, cutoff));
    //         data[i * N_pixel + j] += I * value;
    //     }
    // }
}

void cmos_detection(double input[], double output[], unsigned int data_size)
{
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng * r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    for (unsigned int i(0); i < data_size; ++i)
    {
        const double photons(input[i]);
        const double photoelectrons(
            cmos_detection_function(r, photons));
        // data[i] = photoelectrons;

        output[i] = static_cast<double>(
            cmos_convert_analog_to_digital(photoelectrons));
    }
    gsl_rng_free(r);
}

void emccd_detection(double input[], double output[], unsigned int data_size)
{
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng * r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));

    for (unsigned int i(0); i < data_size; ++i)
    {
        const double photons(input[i]);
        const double photoelectrons(
            emccd_detection_function(r, photons));
        // data[i] = photoelectrons;

        output[i] = static_cast<double>(
            emccd_convert_analog_to_digital(photoelectrons));
    }
    gsl_rng_free(r);
}

inline void detection(double input[], double output[], unsigned int data_size)
{
    cmos_detection(input, output, data_size);
}

} // microscope

#endif /* __MICROSCOPE__MICROSCOPE */

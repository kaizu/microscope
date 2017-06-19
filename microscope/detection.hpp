#ifndef __MICROSCOPE__DETECTION
#define __MICROSCOPE__DETECTION

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>
#include <cmath>

namespace microscope
{

static const double NDist_F40_table[194] = {
    8.313996e-09,
    2.517729e-05,
    1.189114e-03,
    1.830121e-02,
    1.001912e-01,
    2.010323e-01,
    1.918622e-01,
    1.230914e-01,
    7.286084e-02,
    4.673354e-02,
    3.239737e-02,
    2.370227e-02,
    1.798721e-02,
    1.451317e-02,
    1.198214e-02,
    1.037112e-02,
    9.110605e-03,
    8.127194e-03,
    7.294384e-03,
    6.580576e-03,
    6.066170e-03,
    5.498063e-03,
    5.131059e-03,
    4.710254e-03,
    4.393851e-03,
    4.126248e-03,
    3.798044e-03,
    3.580041e-03,
    3.390939e-03,
    3.179437e-03,
    2.963934e-03,
    2.827333e-03,
    2.681631e-03,
    2.525029e-03,
    2.377427e-03,
    2.225326e-03,
    2.089424e-03,
    2.004723e-03,
    1.856721e-03,
    1.828521e-03,
    1.708620e-03,
    1.609919e-03,
    1.503117e-03,
    1.450817e-03,
    1.383616e-03,
    1.328815e-03,
    1.243214e-03,
    1.201214e-03,
    1.138313e-03,
    1.109113e-03,
    9.945515e-04,
    9.854013e-04,
    9.108405e-04,
    8.924303e-04,
    8.371496e-04,
    8.011292e-04,
    7.637088e-04,
    7.379285e-04,
    6.938380e-04,
    6.532975e-04,
    6.097170e-04,
    6.067570e-04,
    5.662965e-04,
    5.424462e-04,
    5.297261e-04,
    4.859556e-04,
    4.597053e-04,
    4.417751e-04,
    4.152348e-04,
    4.248849e-04,
    3.815344e-04,
    3.532441e-04,
    3.484040e-04,
    3.513540e-04,
    3.204037e-04,
    2.911634e-04,
    2.858333e-04,
    2.789232e-04,
    2.608630e-04,
    2.543729e-04,
    2.274426e-04,
    2.273426e-04,
    2.200025e-04,
    2.044424e-04,
    1.961123e-04,
    1.929322e-04,
    1.839821e-04,
    1.778720e-04,
    1.663519e-04,
    1.677119e-04,
    1.516917e-04,
    1.453117e-04,
    1.315915e-04,
    1.273515e-04,
    1.260115e-04,
    1.173414e-04,
    1.104613e-04,
    1.054312e-04,
    1.037812e-04,
    1.009712e-04,
    9.407508e-05,
    8.892702e-05,
    8.411997e-05,
    7.887591e-05,
    7.480286e-05,
    7.094982e-05,
    6.655077e-05,
    6.341173e-05,
    6.146071e-05,
    5.993169e-05,
    5.840867e-05,
    5.651065e-05,
    5.433763e-05,
    5.181960e-05,
    4.911757e-05,
    4.657354e-05,
    4.431651e-05,
    4.224749e-05,
    4.005146e-05,
    3.781244e-05,
    3.571541e-05,
    3.397839e-05,
    3.247937e-05,
    3.109436e-05,
    2.971234e-05,
    2.826733e-05,
    2.698531e-05,
    2.580030e-05,
    2.474528e-05,
    2.390428e-05,
    2.314527e-05,
    2.245626e-05,
    2.177125e-05,
    2.119024e-05,
    2.060724e-05,
    2.001023e-05,
    1.949122e-05,
    1.893922e-05,
    1.836821e-05,
    1.771720e-05,
    1.712720e-05,
    1.646719e-05,
    1.571918e-05,
    1.505817e-05,
    1.434217e-05,
    1.355816e-05,
    1.289915e-05,
    1.221014e-05,
    1.148313e-05,
    1.090313e-05,
    1.031012e-05,
    9.760212e-06,
    9.182206e-06,
    8.768301e-06,
    8.331696e-06,
    7.871391e-06,
    7.548487e-06,
    7.196583e-06,
    6.864779e-06,
    6.495075e-06,
    6.241872e-06,
    5.951469e-06,
    5.626765e-06,
    5.396462e-06,
    5.136059e-06,
    4.885856e-06,
    4.644953e-06,
    4.378450e-06,
    4.198648e-06,
    3.993646e-06,
    3.801944e-06,
    3.590341e-06,
    3.451540e-06,
    3.294638e-06,
    3.149236e-06,
    3.013935e-06,
    2.888633e-06,
    2.749032e-06,
    2.664631e-06,
    2.564630e-06,
    2.471728e-06,
    2.386427e-06,
    2.307127e-06,
    2.232426e-06,
    2.149425e-06,
    2.100924e-06,
    2.042324e-06,
    1.987623e-06,
    1.936822e-06,
    1.889822e-06,
    1.845421e-06,
    1.803921e-06,
    1.764020e-06,
    1.726420e-06
};

double cmos_detection_function(const gsl_rng* rng, const double photons)
{
    const double QE(0.73); // Quantum Efficiency
    const double background_noise(2.0);

    unsigned int signal(
        gsl_ran_poisson(rng, QE * photons + background_noise));

    gsl_ran_discrete_t *f;
    f = gsl_ran_discrete_preproc(194, NDist_F40_table);
    size_t k = gsl_ran_discrete(rng, f);
    const double noise(0.6 + k * 0.1);
    gsl_ran_discrete_free(f);

    const double photoelectrons(
        static_cast<double>(signal) + noise);
    return photoelectrons;
}

double emccd_detection_function(const gsl_rng* rng, const double photons)
{
    const double QE(0.92); // Quantum Efficiency
    const double background_noise(1.0);
    const double EM_gain(300.0);
    const double readout_noise(330.0); //XXX:
    // const double readout_noise(100.0);
    const unsigned int dynodes(11);  // the number of dynode stages

    const double expectation(QE * photons + background_noise);

    double electrons(0.0);

    {
        // Truncated Gaussian Approximation
        // (Eq. 7) in (R. J. Stokey & P. J. Lee, 1983)
        const double alpha(expectation);
        const double A(EM_gain);
        const double nu(static_cast<double>(dynodes));
        const double B(0.5 * (A - 1) / (pow(A, 1 / nu) - 1));
        const double c(exp(alpha * expm1(-A / B)));  // expm1(x) == exp(x) - 1

        if (gsl_rng_uniform(rng) > c)
        {
            const double my(alpha * A);
            const double mx(my / (1 - c));
            const double vary(alpha * A * (A + 2 * B));
            const double varx(vary / (1 - c) - c * mx * mx);

            do
            {
                electrons = mx + gsl_ran_gaussian(rng, sqrt(varx));
            }
            while (electrons < 0);
        }
    }

    const unsigned int signal = std::ceil(electrons);

    // const double alpha(1.0 / EM_gain);
    // double PDF[12000];
    // PDF[0] = exp(-expectation);
    // for (unsigned int electrons = 1; electrons < 12000; ++electrons)
    // {
    //     PDF[electrons] = sqrt(alpha * expectation / electrons) * exp(-alpha * electrons - expectation) * gsl_sf_bessel_I1(2 * sqrt(alpha * expectation * electrons));
    // }

    // gsl_ran_discrete_t *f;
    // f = gsl_ran_discrete_preproc(12000, PDF);
    // const unsigned int signal = gsl_ran_discrete(rng, f);

    // const unsigned int signal = std::ceil(expectation);

    const double noise(gsl_ran_gaussian(rng, readout_noise));
    const double photoelectrons(
        static_cast<double>(signal) + noise);
    // const double photoelectrons(static_cast<double>(signal));
    return photoelectrons;
}

unsigned int convert_analog_to_digital_with_no_fixed_pattern_noise(
    const double photoelectrons,
    const double fullwell, const unsigned int bit, const double offset)
{
    const unsigned int ADC_max(pow(2, bit) - 1);

    // FPN_type == None
    const double gain((fullwell - 0.0) / (ADC_max - offset));
    // const double gain((fullwell - 0.0) / (pow(2.0, bit) - offset));

    const double ADC(
        std::min(photoelectrons, fullwell) / gain + offset);
    return std::max(
        std::min(static_cast<unsigned int>(ADC), ADC_max), (unsigned int)(0));
}

inline unsigned int cmos_convert_analog_to_digital(const double photoelectrons)
{
    const double fullwell(30000);
    const double ADC0(100);
    const unsigned int bit(16);
    const double offset(ADC0);
    return convert_analog_to_digital_with_no_fixed_pattern_noise(photoelectrons, fullwell, bit, offset);
}

inline unsigned int emccd_convert_analog_to_digital(const double photoelectrons)
{
    const double fullwell(370000);
    const unsigned int bit(16);
    const double offset(1743); //XXX
    // const double offset(2000);
    return convert_analog_to_digital_with_no_fixed_pattern_noise(photoelectrons, fullwell, bit, offset);
}


} // microscope

#endif /* __MICROSCOPE__DETECTION */

#ifndef __MICROSCOPE__EMISSION
#define __MICROSCOPE__EMISSION

namespace microscope
{

double tirf_emission(
    double const z,
    double const ATsq = 1.16737230263e+25, double const d = 1.01896194663e+03,
    double const epsilon = 1.0, double const absorption_cross_section = 2e-18,
    double const exposure_time = 100e-3)
{
    // const double ATsq(1.16737230263e+25);
    // const double d(1.01896194663e+03); // nm
    // // const double epsilon(1e-6);
    // const double epsilon(1.0);
    // const double absorption_cross_section(2e-18);
    // // const double absorption_cross_section(1e-9 * 1e-9);
    // const double exposure_time(100e-3);

    return (ATsq * exp(-z / d) * absorption_cross_section
        * exposure_time * epsilon / (4 * M_PI));
}

} // microscope

#endif /* __MICROSCOPE__EMISSION */

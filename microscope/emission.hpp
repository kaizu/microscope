#ifndef __MICROSCOPE__EMISSION
#define __MICROSCOPE__EMISSION

#include <math.h>
#include <utility>

namespace microscope
{

std::pair<double, double> tirf_evanescent_field(
    double const angle, double const wave_length, double const flux_density)
{
    // illumination: assume that uniform illumination (no gaussian)
    // flux density [W/cm2 (joules/sec/cm2)]
    double const P_0 = flux_density * 1e+4;

    // single photon energy
    // hc: (plank const) * (speed of light) [joules meter]
    // wave_length = self.configs.source_wavelength*1e-9
    double const hc = 2.00e-25;
    double const E_wl = hc * 1e+9 / wave_length;

    // photon flux density [photons/sec/cm2]
    double const N_0 = P_0 / E_wl;

    // incident beam amplitude
    double const A2_Is = N_0;
    double const A2_Ip = N_0;

    // incident beam angle
    double const theta = angle / 180 * M_PI;
    double const sin_theta = sin(theta);
    double const cos_theta = sin(theta);
    double const sin2 = sin_theta * sin_theta;
    double const cos2 = cos_theta * cos_theta;

    // index of refraction
    double const n_1 = 1.46;  // fused silica
    double const n_2 = 1.384; // cell
    double const n_3 = 1.337; // culture medium

    double const r = n_2 / n_1;
    double const r2 = r * r;

    assert(sin2 / r2 >= 1);

    double const A2_x = A2_Ip * (4 * cos2 * (sin2 - r2) / (r2 * r2 * cos2 + sin2 - r2));
    double const A2_y = A2_Is * (4 * cos2 / (1 - r2));
    double const A2_z = A2_Ip * (4 * cos2 * sin2 / (r2 * r2 * cos2 + sin2 - r2));
    double const A2_Tp = A2_x + A2_z;
    double const A2_Ts = A2_y;

    double const penetration_depth = wave_length / (4 * M_PI * sqrt(n_1 * n_1 * sin2 - n_2 * n_2));
    double const amplitude = (A2_Tp + A2_Ts) / 2;
    return std::make_pair(amplitude, penetration_depth);
}

double tirf_emission(
    double const z,
    double const ATsq = 1.16737230263e+25, double const d = 1.01896194663e+03,
    double const epsilon = 1.0, double const absorption_cross_section = 2e-18,
    double const exposure_time = 100e-3, double const QY = 0.61)
{
    return (ATsq * exp(-z / d) * QY * absorption_cross_section
        * exposure_time * epsilon / (4 * M_PI));
}

} // microscope

#endif /* __MICROSCOPE__EMISSION */

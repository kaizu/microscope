#ifndef __MICROSCOPE__BORN_WOLF_PSF
#define __MICROSCOPE__BORN_WOLF_PSF

#include "integration.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

namespace microscope
{

struct F_born_wolf_psf_params
{
    double alpha_r;
    double psi;
};

double F_born_wolf_psf_real(double rho, void *params)
{
    struct F_born_wolf_psf_params *p
        = (struct F_born_wolf_psf_params *) params;
    return gsl_sf_bessel_J0(p->alpha_r * rho) * cos(p->psi * rho * rho) * rho;
}

double F_born_wolf_psf_imag(double rho, void *params)
{
    struct F_born_wolf_psf_params *p
        = (struct F_born_wolf_psf_params *) params;
    return gsl_sf_bessel_J0(p->alpha_r * rho) * -sin(p->psi * rho * rho) * rho;
}

double born_wolf_psf(const double r, const double z, const double k, const double N_A)
{
    const double alpha(N_A * k);
    const double alpha_r(alpha * r);
    const double psi(0.5 * alpha * z * N_A);

    struct F_born_wolf_psf_params params = {alpha_r, psi};
    const double result_real(
        integrate1d(&F_born_wolf_psf_real, &params, 0.0, 1.0));
    const double result_imag(
        integrate1d(&F_born_wolf_psf_imag, &params, 0.0, 1.0));
    return result_real * result_real + result_imag * result_imag;
}

}// microscope

#endif /* __MICROSCOPE__BORN_WOLF_PSF */

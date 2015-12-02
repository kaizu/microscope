#ifndef __MICROSCOPE__POINT_SPREADING_FUNCTIONS
#define __MICROSCOPE__POINT_SPREADING_FUNCTIONS

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

#include "integration.hpp"

#include "born_wolf_psf.hpp"
#include "born_wolf_psf_table.hpp"


namespace microscope
{

inline double born_wolf_psf_normalizing_constant(
    const double z, const double k, const double N_A)
{
    // return gsl_pow_2(k * N_A * N_A / (2 * M_PI));
    return gsl_pow_2(k * N_A) / M_PI; //XXX:
}

double born_wolf_psf_simpson(
    const double r, const double z, const double k, const double N_A)
{
    const double alpha(N_A * k);
    const double alpha_r(alpha * r);
    const double psi(0.5 * alpha * z * N_A);

    struct F_born_wolf_psf_params params = {alpha_r, psi};
    const double result_real(
        integrate1d_simpson(&F_born_wolf_psf_real, &params, 0.0, 1.0));
    const double result_imag(
        integrate1d_simpson(&F_born_wolf_psf_imag, &params, 0.0, 1.0));
    return result_real * result_real + result_imag * result_imag;
}

double bilinear_interpolation(
    const double x, const double y,
    const double v00, const double v01, const double v10, const double v11)
{
    return ((1 - x) * ((1 - y) * v00 + y * v01)
        + x * ((1 - y) * v10 + y * v11));
}

double born_wolf_psf_tbl(
    const double r, const double z, const double k, const double N_A)
{
    const double alpha(N_A * k);
    const double alpha_r(alpha * r);
    const double psi(0.5 * alpha * z * N_A);

    const unsigned int N(born_wolf_psf_table::born_wolf_psf_table.N);
    const unsigned int M(born_wolf_psf_table::born_wolf_psf_table.M);
    const double rmax(born_wolf_psf_table::born_wolf_psf_table.rmax);
    const double zmax(born_wolf_psf_table::born_wolf_psf_table.zmax);

    const double n1(alpha_r * N / rmax);
    const double n2(floor(n1));
    const double m1(psi * M / zmax);
    const double m2(floor(m1));
    const double pr(n1 - n2);
    const double pz(m1 - m2);
    const unsigned int n(static_cast<unsigned int>(n2));
    const unsigned int m(static_cast<unsigned int>(m2));

    if (n > N || m > M)
    {
        std::cerr << "Out of bounds." << std::endl;
        return born_wolf_psf(r, z, k, N_A);
    }

    const double
        v00(born_wolf_psf_table::born_wolf_psf_table.y[m * N + n]),
        v10(born_wolf_psf_table::born_wolf_psf_table.y[m * N + n + 1]),
        v01(born_wolf_psf_table::born_wolf_psf_table.y[(m + 1) * N + n]),
        v11(born_wolf_psf_table::born_wolf_psf_table.y[(m + 1) * N + n + 1]);
    return bilinear_interpolation(pr, pz, v00, v01, v10, v11);
}

struct PSF_params
{
    double p[3];
    double k;
    double N_A;
};

double PSF(double const x, double const y, void *params)
{
    struct PSF_params *p = (struct PSF_params *) params;

    const double dx(p->p[0] - x);
    const double dy(p->p[1] - y);
    const double rsq(dx * dx + dy * dy);
    return born_wolf_psf(sqrt(rsq), p->p[2], p->k, p->N_A);
}

double int_psf(
    double const xmin, double const xmax, double const ymin, double const ymax,
    double p[3], double c[3], double const k, double const N_A)
{
    struct PSF_params params = {
        {p[0] - c[0], p[1] - c[1], p[2] - c[2]}, k, N_A};
    const double C(born_wolf_psf_normalizing_constant(p[2] - c[2], k, N_A));
    return C * integrate2d(&PSF, &params, xmin, xmax, ymin, ymax);
}

double PSF_simpson(double const x, double const y, void *params)
{
    struct PSF_params *p = (struct PSF_params *) params;

    const double dx(p->p[0] - x);
    const double dy(p->p[1] - y);
    const double rsq(dx * dx + dy * dy);
    return born_wolf_psf_simpson(sqrt(rsq), p->p[2], p->k, p->N_A);
}

double int_psf_simpson(
    double const xmin, double const xmax, double const ymin, double const ymax,
    double p[3], double c[3], double const k, double const N_A)
{
    struct PSF_params params = {
        {p[0] - c[0], p[1] - c[1], p[2] - c[2]}, k, N_A};
    const double C(born_wolf_psf_normalizing_constant(p[2] - c[2], k, N_A));
    return C * integrate2d_simpson(&PSF, &params, xmin, xmax, ymin, ymax);
    // return C * integrate2d_simpson(&PSF_simpson, &params, xmin, xmax, ymin, ymax);
}

double PSF_tbl(double const x, double const y, void *params)
{
    struct PSF_params *p = (struct PSF_params *) params;

    const double dx(p->p[0] - x);
    const double dy(p->p[1] - y);
    const double rsq(dx * dx + dy * dy);
    return born_wolf_psf_tbl(sqrt(rsq), p->p[2], p->k, p->N_A);
}

double int_psf_tbl(
    double const xmin, double const xmax, double const ymin, double const ymax,
    double p[3], double c[3], double const k, double const N_A)
{
    struct PSF_params params = {
        {p[0] - c[0], p[1] - c[1], p[2] - c[2]}, k, N_A};
    const double C(born_wolf_psf_normalizing_constant(p[2] - c[2], k, N_A));
    return C * integrate2d_simpson(&PSF_tbl, &params, xmin, xmax, ymin, ymax);
}

struct PSF_cylinder_params
{
    double z;
    double k;
    double N_A;
    double v;
};

double PSF_cylinder(double const r, void *params)
{
    struct PSF_cylinder_params *p = (struct PSF_cylinder_params *) params;
    return r * born_wolf_psf(r, p->z, p->k, p->N_A);
}

double int_psf_cylinder(
    double const rmin, double const rmax, PSF_cylinder_params * params)
{
    const double C(born_wolf_psf_normalizing_constant(params->z, params->k, params->N_A));
    return (2 * M_PI) * C * integrate1d(&PSF_cylinder, params, rmin, rmax);
}

double int_psf_cylinder(double const r, PSF_cylinder_params * params)
{
    const double C(born_wolf_psf_normalizing_constant(params->z, params->k, params->N_A));
    return (2 * M_PI) * C * integrate1d(&PSF_cylinder, params, 0.0, r) - params->v;
}

double int_psf_cylinder(double const r, void * params)
{
    return int_psf_cylinder(r, static_cast<PSF_cylinder_params*>(params));
}

double int_psf_cylinder(
    double const rmin, double const rmax, double const z, double const k, double const N_A)
{
    struct PSF_cylinder_params params = {z, k, N_A};
    return int_psf_cylinder(rmin, rmax, &params);
}

double int_psf_gaussian(
    double const xmin, double const xmax, double const ymin, double const ymax,
    double p[3], double c[3], double const k, double const N_A)
{
    // const double sigma(sqrt(2.0) / (k * N_A)); //  The variance of Gaussian changes linearly with the axial axis, p[2].
    const double sigma(
        1.0 / (sqrt(2 * born_wolf_psf_tbl(0.0, p[2] - c[2], k, N_A)) * k * N_A));
    const double sigma_inv(1.0 / (sigma * sqrt(2.0)));
    const double term1(
        gsl_sf_erf((xmax - p[0] + c[0]) * sigma_inv)
        - gsl_sf_erf((xmin - p[0] + c[0]) * sigma_inv));
    const double term2(
        gsl_sf_erf((ymax - p[1] + c[1]) * sigma_inv)
        - gsl_sf_erf((ymin - p[1] + c[1]) * sigma_inv));
    return 0.25 * term1 * term2;
}

} // microscope

#endif /* __MICROSCOPE__POINT_SPREADING_FUNCTIONS */

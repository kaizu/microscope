#ifndef __MICROSCOPE__POINT_SPREADING_FUNCTIONS
#define __MICROSCOPE__POINT_SPREADING_FUNCTIONS

#include <gsl/gsl_math.h>

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

double PSF(double x, double y, void *params)
{
    struct PSF_params *p = (struct PSF_params *) params;

    const double dx(p->p[0] - x);
    const double dy(p->p[1] - y);
    const double rsq(dx * dx + dy * dy);
    return born_wolf_psf(sqrt(rsq), p->p[2], p->k, p->N_A);
}

double int_psf(
    double xmin, double xmax, double ymin, double ymax,
    double p[3], double c[3], double k, double N_A)
{
    struct PSF_params params = {
        {p[0] - c[0], p[1] - c[1], p[2] - c[2]}, k, N_A};
    const double C(gsl_pow_2(k * N_A * N_A / (2 * M_PI)));
    return C * integrate2d(&PSF, &params, xmin, xmax, ymin, ymax);
}

double PSF_simpson(double x, double y, void *params)
{
    struct PSF_params *p = (struct PSF_params *) params;

    const double dx(p->p[0] - x);
    const double dy(p->p[1] - y);
    const double rsq(dx * dx + dy * dy);
    return born_wolf_psf_simpson(sqrt(rsq), p->p[2], p->k, p->N_A);
}

double int_psf_simpson(
    double xmin, double xmax, double ymin, double ymax,
    double p[3], double c[3], double k, double N_A)
{
    struct PSF_params params = {
        {p[0] - c[0], p[1] - c[1], p[2] - c[2]}, k, N_A};
    const double C(gsl_pow_2(k * N_A * N_A / (2 * M_PI)));
    return C * integrate2d_simpson(&PSF, &params, xmin, xmax, ymin, ymax);
    // return C * integrate2d_simpson(&PSF_simpson, &params, xmin, xmax, ymin, ymax);
}

double PSF_tbl(double x, double y, void *params)
{
    struct PSF_params *p = (struct PSF_params *) params;

    const double dx(p->p[0] - x);
    const double dy(p->p[1] - y);
    const double rsq(dx * dx + dy * dy);
    return born_wolf_psf_tbl(sqrt(rsq), p->p[2], p->k, p->N_A);
}

double int_psf_tbl(
    double xmin, double xmax, double ymin, double ymax,
    double p[3], double c[3], double k, double N_A)
{
    struct PSF_params params = {
        {p[0] - c[0], p[1] - c[1], p[2] - c[2]}, k, N_A};
    const double C(gsl_pow_2(k * N_A * N_A / (2 * M_PI)));
    return C * integrate2d_simpson(&PSF_tbl, &params, xmin, xmax, ymin, ymax);
}

struct PSF_cylinder_params
{
    double z;
    double k;
    double N_A;
};

double PSF_cylinder(double r, void *params)
{
    struct PSF_cylinder_params *p = (struct PSF_cylinder_params *) params;
    return r * born_wolf_psf(r, p->z, p->k, p->N_A);
}

double int_psf_cylinder(
    double rmin, double rmax, double z, double k, double N_A)
{
    struct PSF_cylinder_params params = {z, k, N_A};
    return gsl_pow_2(k * N_A * N_A) * integrate1d(&PSF_cylinder, &params, rmin, rmax);
}

} // microscope

#endif /* __MICROSCOPE__POINT_SPREADING_FUNCTIONS */

#include <iostream>
#include <gsl/gsl_rng.h>
#include <boost/array.hpp>

#include "point_spreading_functions.hpp"
using namespace microscope;


void test_point_spreading_function_cutoff(double k, double N_A)
{
    const double alpha_inv(1.0 / (k * N_A));
    const double psi_inv(2 * alpha_inv / N_A);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(16);

    std::cout << "alpha = " << 1.0 / alpha_inv << std::endl;
    std::cout << "1 um corresponds to " << 1000 / alpha_inv << " / alpha" << std::endl;
    const double z(0.0 * psi_inv * 0.5);
    const double I0(int_psf_cylinder(0, 1000, z, k, N_A));
    const double I1(int_psf_cylinder(0, 10000, z, k, N_A));
    const double I2(int_psf_cylinder(0, 30000, z, k, N_A));
    const double I3(int_psf_cylinder(0, 50000, z, k, N_A));
    std::cout << "The integration over r < 1 um = " << I0 << std::endl;
    std::cout << "An integration over r < 10 um = " << I1 << std::endl;
    std::cout << "An integration over r < 30 um = " << I2 << std::endl;
    std::cout << "An integration over r < 50 um = " << I3 << std::endl;
    std::cout << "The relative error  1-50 um = " << (I3 - I0) / I3 << std::endl;
    std::cout << "The relative error 10-50 um = " << (I3 - I1) / I3 << std::endl;
    std::cout << "The relative error 30-50 um = " << (I3 - I2) / I3 << std::endl;
}

double int_psf_with_cutoff(
    double p[3], double c[3], double k, double N_A,
    double xmin, double xmax, double ymin, double ymax,
    double cutoff)
{
    const double rsq_min(
        std::min(gsl_pow_2(p[0] - c[0] - xmin), gsl_pow_2(p[0] - c[0] - xmax))
        + std::min(gsl_pow_2(p[1] - c[1] - ymin), gsl_pow_2(p[1] - c[1] - ymax)));
    if (rsq_min >= cutoff * cutoff)
    {
        return 0.0;
    }
    else
    {
        const double result(int_psf(xmin, xmax, ymin, ymax, p, c, k, N_A));
        return result;
    }
}

double overlay_psf(
    double data[], unsigned int N_pixel, double pixel_length,
    double p[3], double c[3], double k, double N_A,
    double cutoff)
{
    for (unsigned int i(0); i < N_pixel; ++i)
    {
        const double xmin(pixel_length * (i - N_pixel * 0.5));
        const double xmax(xmin + pixel_length);
        for (unsigned int j(0); j < N_pixel; ++j)
        {
            const double ymin(pixel_length * (j - N_pixel * 0.5));
            const double ymax(ymin + pixel_length);

            const double value(int_psf_with_cutoff(
                p, c, k, N_A, xmin, xmax, ymin, ymax, cutoff));
            data[i * N_pixel + j] += value;
        }
    }
}

void generate_random_points(
    double points[][3], const unsigned int N_point, const double L)
{
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng * r = gsl_rng_alloc(T);

    for (unsigned int i(0); i < N_point; ++i)
    {
        points[i][0] = L * (0.5 - gsl_rng_uniform(r));
        points[i][1] = L * (0.5 - gsl_rng_uniform(r));
        points[i][2] = 0.0;
    }

    gsl_rng_free(r);
}

int main()
{
    const double lambda(508); // nm
    const double N_A(1.4);
    const double k(2 * M_PI / lambda);

    // test_point_spreading_function_cutoff(k, N_A);

    const unsigned int N_pixel(600);
    const double pixel_length(6500 / 100);
    const double L(pixel_length * N_pixel);
    const double L_2(L * 0.5);

    // double p[] = {0, 0, 0};
    double focal_point[] = {0, 0, 0};

    const unsigned int N_point(10);
    double points[N_point][3];
    generate_random_points(points, N_point, L);

    boost::array<double, N_pixel * N_pixel> data;
    data.fill(0.0);

    const double cutoff(1000);
    for (unsigned int i(0); i < N_point; ++i)
    {
        overlay_psf(data.data(), N_pixel, pixel_length,
            points[i], focal_point, k, N_A, cutoff);
    }

    std::cout.setf(std::ios::scientific);
    std::cout.precision(16);
    for (unsigned int i(0); i < N_pixel; ++i)
    {
        for (unsigned int j(0); j < N_pixel; ++j)
        {
            std::cout << data[i * N_pixel + j] << std::endl;
        }
    }

    // for (unsigned int i(0); i < 200; ++i)
    // {
    //     const double r(20 * alpha_inv * i / 200);
    //     for (unsigned int j(0); j < 200; ++j)
    //     {
    //         const double z(20 * alpha_inv * j / 200);
    //         std::cout << r << "," << z << ","
    //             << born_wolf_psf(r, z, k, N_A) << std::endl;
    //     }
    // }

    // double p[] = {0, 0, 0};
    // double c[] = {0, 0, 0};
    // const unsigned int N_pixel(20);
    // // const unsigned int N_pixel(600);
    // // const double pixel_length(65);
    // const double pixel_length(40 * alpha_inv / N_pixel);
    // const double L(N_pixel * pixel_length);

    // double tot(0.0);
    // for (unsigned int i(0); i < N_pixel; ++i)
    // {
    //     const double xmin(pixel_length * (i - N_pixel * 0.5));
    //     const double xmax(xmin + pixel_length);
    //     for (unsigned int j(0); j < N_pixel; ++j)
    //     {
    //         const double ymin(pixel_length * (j - N_pixel * 0.5));
    //         const double ymax(ymin + pixel_length);

    //         const double result(int_psf(xmin, xmax, ymin, ymax, p, c, k, N_A));
    //         tot += result;
    //         std::cout << i << "," << j << "," << result << std::endl;
    //     }
    // }

    // std::cout.precision(16);
    // std::cout << tot << std::endl;
    // std::cout << int_psf(-0.5 * L, +0.5 * L, -0.5 * L, +0.5 * L, p, c, k, N_A) << std::endl;

    // double cutoff;
    // std::cout << int_psf_cylinder(
    //     0, alpha_inv * cutoff, -cutoff * psi_inv, +cutoff * psi_inv, k, N_A) << std::endl;
    // std::cout.precision(16);
    // double prev(0);
    // for (unsigned int i(1); i < 1000; ++i)
    // {
    //     cutoff = static_cast<double>(i);
    //     struct func_params params = {k, N_A, cutoff * psi_inv};
    //     const double val(gsl_pow_2(k * N_A * N_A) * integrate1d(&func, &params, 0, 10 * cutoff * alpha_inv));
    //     // const double val(gsl_pow_2(k * N_A * N_A) * integrate1d(&func, &params, 0, cutoff * alpha_inv));
    //     std::cout << cutoff * psi_inv << "," << val << "," << (val - prev) / val << std::endl;
    //     prev = val;
    //     // std::cout << i << "," << int_psf_cylinder_new(
    //     //     0, alpha_inv * cutoff, -20 * psi_inv, +20 * psi_inv, k, N_A) << std::endl;
    // }

    //cutoff = 200;
    //std::cout << int_psf_cylinder_new(
    //    0, alpha_inv * cutoff, -cutoff * psi_inv, +cutoff * psi_inv, k, N_A) << std::endl;
    // cutoff = 90;
    // std::cout << int_psf_cylinder(
    //     0, alpha_inv * cutoff, -cutoff * psi_inv, +cutoff * psi_inv, k, N_A) << std::endl;
    // cutoff = 100;
    // std::cout << int_psf_cylinder(
    //     0, alpha_inv * cutoff, -cutoff * psi_inv, +cutoff * psi_inv, k, N_A) << std::endl;
    // cutoff = 120;
    // std::cout << int_psf_cylinder(
    //     0, alpha_inv * cutoff, -cutoff * psi_inv, +cutoff * psi_inv, k, N_A) << std::endl;
    // std::cout << int_psf_cylinder(
    //     0, alpha_inv * 40, -40 * psi_inv, +40 * psi_inv, k, N_A) << std::endl;

    // struct func_params params = {k, N_A, 0.0};
    // std::cout << integrate1d(&func, &params, 0, 999) << std::endl;
    // std::cout << M_PI * (alpha_inv * 40 * alpha_inv * 40) * 80 * psi_inv << std::endl;
}

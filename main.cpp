#include <iostream>

#include "point_spreading_functions.hpp"
using namespace microscope;


int main()
{
    const double lambda(508); // nm
    const double N_A(1.4);
    const double k(2 * M_PI / lambda);
    const double alpha_inv(1.0 / (k * N_A));
    const double psi_inv(2 * alpha_inv / N_A);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(16);

    std::cout << int_psf_cylinder(0, 10 * alpha_inv, 0, k, N_A) << std::endl;

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

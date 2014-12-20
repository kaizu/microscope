#include <iostream>

#include "point_spreading_functions.hpp"
using namespace microscope;


int main()
{
    const double lambda(508); // nm
    const double N_A(1.4);
    const double k(2 * M_PI / lambda);
    const double alpha_inv(lambda / (2 * M_PI * N_A));

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

    double p[] = {0, 0, 0};
    double c[] = {0, 0, 0};
    const unsigned int N_pixel(20);
    // const unsigned int N_pixel(600);
    // const double pixel_length(65);
    const double pixel_length(40 * alpha_inv / N_pixel);
    const double L(N_pixel * pixel_length);

    double tot(0.0);
    for (unsigned int i(0); i < N_pixel; ++i)
    {
        const double xmin(pixel_length * (i - N_pixel * 0.5));
        const double xmax(xmin + pixel_length);
        for (unsigned int j(0); j < N_pixel; ++j)
        {
            const double ymin(pixel_length * (j - N_pixel * 0.5));
            const double ymax(ymin + pixel_length);

            const double result(int_psf(xmin, xmax, ymin, ymax, p, c, k, N_A));
            tot += result;
            std::cout << i << "," << j << "," << result << std::endl;
        }
    }

    std::cout.precision(16);
    std::cout << tot << std::endl;
    std::cout << int_psf(-0.5 * L, +0.5 * L, -0.5 * L, +0.5 * L, p, c, k, N_A) << std::endl;
}

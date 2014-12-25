#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_rng.h>

#include <boost/array.hpp>
#include <boost/range/numeric.hpp>

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
        const double result(int_psf_simpson(xmin, xmax, ymin, ymax, p, c, k, N_A));
        // const double result(int_psf(xmin, xmax, ymin, ymax, p, c, k, N_A));
        return result;
    }
}

double emission(double z)
{
    const double ATsq(1.16737230263e+25);
    const double d(1.01896194663e+03); // nm
    const double epsilon(1e-6);
    const double cross_section(1e-9 * 1e-9);
    const double exposure_time(100e-3);
    return ATsq * exp(-z / d) * cross_section * exposure_time * epsilon / (4 * M_PI);
}

void overlay_psf(
    double data[], unsigned int N_pixel, double pixel_length,
    double p[3], double I, double c[3], double k, double N_A,
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
            data[i * N_pixel + j] += I * value;
        }
    }
}

void generate_random_points(
    double points[][3], double intensity[], const unsigned int N_point, const double L)
{
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng * r = gsl_rng_alloc(T);

    for (unsigned int i(0); i < N_point; ++i)
    {
        points[i][0] = L * (0.5 - gsl_rng_uniform(r));
        points[i][1] = L * (0.5 - gsl_rng_uniform(r));
        points[i][2] = gsl_rng_uniform(r) * 400;
        // points[i][2] = 0.0;
        intensity[i] = emission(points[i][2]);
    }

    gsl_rng_free(r);
}

void save_data(const char filename[], double data[], unsigned int data_size)
{
    std::ofstream fout;
    fout.open(filename);
    fout.setf(std::ios::scientific);
    fout.precision(16);
    for (unsigned int i(0); i < data_size; ++i)
    {
        fout << data[i] << std::endl;
    }
    fout.close();
}

std::string point_as_str(double p[3])
{
    std::stringstream sout;
    sout << std::showpos;
    sout << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return sout.str();
}

int main()
{
    const double lambda(508); // nm
    const double N_A(1.4);
    const double k(2 * M_PI / lambda);

    // test_point_spreading_function_cutoff(k, N_A);

    const unsigned int N_pixel(60);
    const double pixel_length(6500 / 100);
    const double L(pixel_length * N_pixel);
    const double L_2(L * 0.5);
    double focal_point[] = {0, 0, 0};

    const unsigned int N_point(10);
    double points[N_point][3];
    double intensity[N_point];
    generate_random_points(points, intensity, N_point, L);

    boost::array<double, N_pixel * N_pixel> data;
    data.fill(0.0);

    const double cutoff(1000);
    double Itot(0.0);
    for (unsigned int i(0); i < N_point; ++i)
    {
        std::cout << "[" << i << "] = " << point_as_str(points[i]) << std::endl;

        overlay_psf(data.data(), N_pixel, pixel_length,
            points[i], intensity[i], focal_point, k, N_A, cutoff);

        Itot += intensity[i] * int_psf(
            -L_2, +L_2, -L_2, +L_2, points[i], focal_point, k, N_A);
    }

    save_data("result.txt", data.data(), N_pixel * N_pixel);

    const double I(boost::accumulate(data, 0.0));
    std::cout << "The total intensity of an output image is " << I << std::endl;
    std::cout << "The total intensity of PSFs is " << Itot << std::endl;
    std::cout << "The relative error is " << fabs(Itot - I) / Itot << std::endl;
}

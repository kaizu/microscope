#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_rng.h>

#include <boost/array.hpp>
#include <boost/range/numeric.hpp>
#include <boost/format.hpp>

#include "point_spreading_functions.hpp"
#include "detection.hpp"
using namespace microscope;


std::string point_as_str(double p[3])
{
    std::stringstream sout;
    sout << std::showpos;
    sout << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return sout.str();
}

double emission(double z)
{
    const double ATsq(1.16737230263e+25);
    const double d(1.01896194663e+03); // nm
    // const double epsilon(1e-6);
    const double epsilon(1.0);
    const double absorption_cross_section(2e-18);
    // const double absorption_cross_section(1e-9 * 1e-9);
    const double exposure_time(100e-3);

    return (ATsq * exp(-z / d) * absorption_cross_section
        * exposure_time * epsilon / (4 * M_PI));
}

void overlay_psf(
    double data[], unsigned int N_pixel, double pixel_length,
    double p[3], double I, double c[3], double k, double N_A,
    double cutoff)
{
    const double offset(N_pixel * pixel_length * -0.5);
    const double x(p[0] - c[0] - offset);
    const double y(p[1] - c[1] - offset);
    const unsigned int imin(
        static_cast<unsigned int>(
            std::max(0.0, floor((x - cutoff) / pixel_length))));
    const unsigned int imax(
        std::min(N_pixel,
            static_cast<unsigned int>(ceil((x + cutoff) / pixel_length))));
    const unsigned int jmin(
        static_cast<unsigned int>(
            std::max(0.0, floor((y - cutoff) / pixel_length))));
    const unsigned int jmax(
        std::min(N_pixel,
            static_cast<unsigned int>(ceil((y + cutoff) / pixel_length))));

    for (unsigned int i(imin); i < imax; ++i)
    {
        const double xmin(pixel_length * i + offset);
        const double xmax(xmin + pixel_length);

        for (unsigned int j(jmin); j < jmax; ++j)
        {
            const double ymin(pixel_length * j + offset);
            const double ymax(ymin + pixel_length);

            const double value(int_psf_tbl(
                xmin, xmax, ymin, ymax, p, c, k, N_A));
            data[i * N_pixel + j] += I * value;
        }
    }

    // for (unsigned int i(0); i < N_pixel; ++i)
    // {
    //     const double xmin(pixel_length * (i - N_pixel * 0.5));
    //     const double xmax(xmin + pixel_length);
    //     for (unsigned int j(0); j < N_pixel; ++j)
    //     {
    //         const double ymin(pixel_length * (j - N_pixel * 0.5));
    //         const double ymax(ymin + pixel_length);

    //         const double value(int_psf_with_cutoff(
    //             p, c, k, N_A, xmin, xmax, ymin, ymax, cutoff));
    //         data[i * N_pixel + j] += I * value;
    //     }
    // }
}

void generate_random_points(
    double points[][3], double intensity[], const unsigned int N_point, const double L)
{
    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng * r = gsl_rng_alloc(T);

    for (unsigned int i(0); i < N_point - 200; ++i)
    {
        points[i][0] = L * (0.5 - gsl_rng_uniform(r));
        points[i][1] = L * (0.5 - gsl_rng_uniform(r));
        points[i][2] = gsl_rng_uniform(r) * 3000;
        // points[i][2] = 0.0;
        intensity[i] = emission(points[i][2]);
    }

    for (unsigned int i(N_point - 200); i < N_point; ++i)
    {
        points[i][0] = L * (0.5 - gsl_rng_uniform(r));
        points[i][1] = L * (0.5 - gsl_rng_uniform(r));
        points[i][2] = 0.0;
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

void read_input(const char filename[], double points[][3], double intensity[],
                unsigned int data_size)
{
    std::ifstream fin(filename);
    std::string buf;

    std::getline(fin, buf); // Ignore a header line

    unsigned int i(0);
    while (fin && std::getline(fin, buf) && i < data_size)
    {
        std::string token;
        std::istringstream stream(buf);
        std::stringstream ss;
        double val;
        unsigned int sid;

        double x, y, z;

        {
            std::getline(stream, token, ',');
            ss << token;
            ss >> val;
            x = val;
            ss.str("");
            ss.clear(std::stringstream::goodbit);
        }
        {
            std::getline(stream, token, ',');
            ss << token;
            ss >> val;
            y = val;
            ss.str("");
            ss.clear(std::stringstream::goodbit);
        }
        {
            std::getline(stream, token, ',');
            ss << token;
            ss >> val;
            z = val;
            ss.str("");
            ss.clear(std::stringstream::goodbit);
        }
        {
            std::getline(stream, token, ',');
            std::getline(stream, token, ',');
            ss << token;
            ss >> sid;
            ss.str("");
            ss.clear(std::stringstream::goodbit);
        }

        // if (sid != 0)
        // {
        //     continue;
        // }

        points[i][0] = (x - 15) * 1000;
        points[i][1] = (y - 15) * 1000;
        points[i][2] = fabs((z - 1.5) * 1000);
        // points[i][2] = z * 1000;
        intensity[i] = emission(points[i][2]);
        // std::cout << point_as_str(points[i]) << std::endl;
        i++;
    };
}

int main()
{
    const double lambda(508); // nm
    const double N_A(1.4);
    const double k(2 * M_PI / lambda);

    const unsigned int N_pixel(600);
    const double pixel_length(6500 / 100);
    const double L(pixel_length * N_pixel);
    const double L_2(L * 0.5);
    double focal_point[] = {0, 0, 0};

    // const unsigned int N_point(17200);
    // const unsigned int N_point(1000);
    const unsigned int N_point(2620);
    double points[N_point][3];
    double intensity[N_point];

    boost::array<double, N_pixel * N_pixel> data;
    data.fill(0.0);

    // const double cutoff(20000);
    const double cutoff(2000);
    double Itot(0.0);

    if (false)
    {
        for (unsigned int i = 0; i < 1001; ++i)
        {
            double const z = static_cast<double>(i) / (200 * k * N_A);
            double const val = born_wolf_psf(0.0, z, k, N_A);
            std::cout << z << "\t" << 1.0 / gsl_pow_2(4 * val) << std::endl;
            // std::cout << z << "\t" << 0.25 / (val * val) / (k * N_A) << std::endl;
        }
    }

    // if (false)
    {
        double const z = 4.0 / (k * N_A);
        double const C = born_wolf_psf_normalizing_constant(z, k, N_A);
        for (unsigned int i = 0; i < 2001; ++i)
        {
            double const r = (static_cast<double>(i) - 1000) / (150 * k * N_A);
            std::cout << r
                << "\t" << C * born_wolf_psf(r, 0.0, k, N_A)
                << "\t" << C * born_wolf_psf(r, z, k, N_A) << std::endl;
        }
    }

    if (false)
    {
        gsl_rng_env_setup();
        const gsl_rng_type * T = gsl_rng_default;
        gsl_rng * r = gsl_rng_alloc(T);

        const gsl_root_fsolver_type* solverType(gsl_root_fsolver_brent);
        gsl_root_fsolver* solver(gsl_root_fsolver_alloc(solverType));

        struct PSF_cylinder_params params = {1000.0, k, N_A, 0.0};
        const double vmax = int_psf_cylinder(4000, &params);
        double rv(0.0);
        for (unsigned int i = 0; i < 101; ++i)
        {
            params.v = vmax * i / 100;
            rv = find_root(&int_psf_cylinder, &params, solver, rv, 4000, 1e-18, 1e-12);
            std::cout << i << "\t" << rv << "\t"
                << int_psf_cylinder(rv, &params) + params.v << std::endl;
        }

        gsl_root_fsolver_free(solver);
        gsl_rng_free(r);

    }

    // generate_random_points(points, intensity, N_point, L);

    // for (unsigned int i(0); i < N_point; ++i)
    // {
    //     std::cout << "[" << i << "] = " << point_as_str(points[i]) << std::endl;

    //     overlay_psf(data.data(), N_pixel, pixel_length,
    //         points[i], intensity[i], focal_point, k, N_A, cutoff);
    //     Itot += intensity[i];
    // }

    // // const unsigned int frames(101);
    // const unsigned int frames(1);
    // for (unsigned int cnt(0); cnt < frames; ++cnt)
    // {
    //     // generate_random_points(points, intensity, N_point, L);
    //     read_input(
    //         (boost::format("sample_data3/test%03d.csv")
    //             % (cnt * (101 / frames))).str().c_str(),
    //         points, intensity, N_point);

    //     for (unsigned int i(0); i < N_point; ++i)
    //     {
    //         std::cout << "[" << cnt << ":" << i << "/"
    //             << N_point << "] = " << point_as_str(points[i]) << std::endl;

    //         const double I(intensity[i] / frames);
    //         // overlay_psf(data.data(), N_pixel, pixel_length,
    //         //     points[i], I, focal_point, k, N_A, cutoff);

    //         const double QE(0.73);
    //         unsigned int signal(gsl_ran_poisson(r, QE * I));
    //         std::cout << signal << std::endl;

    //         Itot += I;
    //     }
    // }

    // const double I(boost::accumulate(data, 0.0));
    // std::cout << "The total intensity of an output image is " << I << std::endl;
    // std::cout << "The total intensity of PSFs is " << Itot << std::endl;
    // std::cout << "The relative error is " << fabs(Itot - I) / Itot << std::endl;

    // {
    //     gsl_rng_env_setup();
    //     const gsl_rng_type * T = gsl_rng_default;
    //     gsl_rng * r = gsl_rng_alloc(T);

    //     for (unsigned int i(0); i < N_pixel; ++i)
    //     {
    //         for (unsigned int j(0); j < N_pixel; ++j)
    //         {
    //             const double photons(data[i * N_pixel + j]);
    //             const double photoelectrons(
    //                 cmos_detection_function(r, photons));
    //             // data[i * N_pixel + j] = photoelectrons;
    //             data[i * N_pixel + j] = static_cast<double>(
    //                 convert_analog_to_digital(photoelectrons));
    //         }
    //     }
    //     gsl_rng_free(r);
    // }

    // save_data((boost::format("result%03d.txt") % (frames)).str().c_str(),
    //           data.data(), N_pixel * N_pixel);
}

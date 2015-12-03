#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_rng.h>

#include <boost/array.hpp>
#include <boost/range/numeric.hpp>
#include <boost/format.hpp>

#include "microscope.hpp"

using namespace microscope;


std::string point_as_str(double p[3])
{
    std::stringstream sout;
    sout << std::showpos;
    sout << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return sout.str();
}

// void test_point_spreading_function_cutoff(double const k, double const N_A)
// {
//     const double alpha_inv(1.0 / (k * N_A));
//     const double psi_inv(2 * alpha_inv / N_A);
// 
//     std::cout.setf(std::ios::scientific);
//     std::cout.precision(16);
// 
//     std::cout << "alpha = " << 1.0 / alpha_inv << std::endl;
//     std::cout << "1 um corresponds to " << 1000 / alpha_inv << " / alpha" << std::endl;
//     const double z(0.0 * psi_inv * 0.5);
//     const double I0(int_psf_cylinder(0, 1000, z, k, N_A));
//     const double I1(int_psf_cylinder(0, 10000, z, k, N_A));
//     const double I2(int_psf_cylinder(0, 30000, z, k, N_A));
//     const double I3(int_psf_cylinder(0, 50000, z, k, N_A));
//     std::cout << "The integration over r < 1 um = " << I0 << std::endl;
//     std::cout << "An integration over r < 10 um = " << I1 << std::endl;
//     std::cout << "An integration over r < 30 um = " << I2 << std::endl;
//     std::cout << "An integration over r < 50 um = " << I3 << std::endl;
//     std::cout << "The relative error  1-50 um = " << (I3 - I0) / I3 << std::endl;
//     std::cout << "The relative error 10-50 um = " << (I3 - I1) / I3 << std::endl;
//     std::cout << "The relative error 30-50 um = " << (I3 - I2) / I3 << std::endl;
// }

// double int_psf_with_cutoff(
//     double p[3], double c[3], double const k, double const N_A,
//     double const xmin, double const xmax, double const ymin, double const ymax,
//     double const cutoff)
// {
//     const double rsq_min(
//         std::min(gsl_pow_2(p[0] - c[0] - xmin), gsl_pow_2(p[0] - c[0] - xmax))
//         + std::min(gsl_pow_2(p[1] - c[1] - ymin), gsl_pow_2(p[1] - c[1] - ymax)));
//     if (rsq_min >= cutoff * cutoff)
//     {
//         return 0.0;
//     }
//     else
//     {
//         const double result(int_psf_tbl(xmin, xmax, ymin, ymax, p, c, k, N_A));
//         // const double result(int_psf_simpson(xmin, xmax, ymin, ymax, p, c, k, N_A));
//         // const double result(int_psf(xmin, xmax, ymin, ymax, p, c, k, N_A));
//         return result;
//     }
// }

// void overlay_psf(
//     double data[], unsigned int const N_pixel, double const pixel_length,
//     double p[3], double const I, double c[3], double const k, double const N_A,
//     double const cutoff)
// {
//     const double offset(N_pixel * pixel_length * -0.5);
//     const double x(p[0] - c[0] - offset);
//     const double y(p[1] - c[1] - offset);
//     const unsigned int imin(
//         static_cast<unsigned int>(
//             std::max(0.0, floor((x - cutoff) / pixel_length))));
//     const unsigned int imax(
//         std::min(N_pixel,
//             static_cast<unsigned int>(ceil((x + cutoff) / pixel_length))));
//     const unsigned int jmin(
//         static_cast<unsigned int>(
//             std::max(0.0, floor((y - cutoff) / pixel_length))));
//     const unsigned int jmax(
//         std::min(N_pixel,
//             static_cast<unsigned int>(ceil((y + cutoff) / pixel_length))));
// 
//     for (unsigned int i(imin); i < imax; ++i)
//     {
//         const double xmin(pixel_length * i + offset);
//         const double xmax(xmin + pixel_length);
// 
//         for (unsigned int j(jmin); j < jmax; ++j)
//         {
//             const double ymin(pixel_length * j + offset);
//             const double ymax(ymin + pixel_length);
// 
//             // const double value(int_psf_gaussian(
//             //     xmin, xmax, ymin, ymax, p, c, k, N_A));
//             const double value(int_psf_tbl(
//                 xmin, xmax, ymin, ymax, p, c, k, N_A));
//             data[i * N_pixel + j] += I * value;
//         }
//     }
// 
//     // for (unsigned int i(0); i < N_pixel; ++i)
//     // {
//     //     const double xmin(pixel_length * (i - N_pixel * 0.5));
//     //     const double xmax(xmin + pixel_length);
//     //     for (unsigned int j(0); j < N_pixel; ++j)
//     //     {
//     //         const double ymin(pixel_length * (j - N_pixel * 0.5));
//     //         const double ymax(ymin + pixel_length);
// 
//     //         const double value(int_psf_with_cutoff(
//     //             p, c, k, N_A, xmin, xmax, ymin, ymax, cutoff));
//     //         data[i * N_pixel + j] += I * value;
//     //     }
//     // }
// }

void generate_random_points(
    double points[][3], unsigned int const N_point, double const L)
{
    assert(N_point >= 200);

    gsl_rng_env_setup();
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng * r = gsl_rng_alloc(T);

    for (unsigned int i(0); i < N_point - 200; ++i)
    {
        points[i][0] = L * (0.5 - gsl_rng_uniform(r));
        points[i][1] = L * (0.5 - gsl_rng_uniform(r));
        points[i][2] = gsl_rng_uniform(r) * 3000;
        // intensity[i] = emission(points[i][2]);
    }

    for (unsigned int i(N_point - 200); i < N_point; ++i)
    {
        points[i][0] = L * (0.5 - gsl_rng_uniform(r));
        points[i][1] = L * (0.5 - gsl_rng_uniform(r));
        points[i][2] = 0.0;
        // intensity[i] = emission(points[i][2]);
    }

    gsl_rng_free(r);
}

void save_data(char const filename[], double data[], unsigned int const data_size)
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

void read_input(
    char const filename[], double points[][3], unsigned int const data_size,
    double shift[3], double const scale)
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

        points[i][0] = (x - shift[0]) * scale;
        points[i][1] = (y - shift[1]) * scale;
        points[i][2] = fabs(z - shift[2]) * scale;
        // intensity[i] = emission(points[i][2]);
        // std::cout << point_as_str(points[i]) << std::endl;
        i++;
    };
}

int main(int argc, char** argv)
{
    const double lambda(508); // nm
    const double N_A(1.4);
    const double k(2 * M_PI / lambda);

    // test_point_spreading_function_cutoff(k, N_A);

    const unsigned int N_pixel(600);
    const double pixel_length(6500 / 100);
    const double L(pixel_length * N_pixel);
    const double L_2(L * 0.5);
    double focal_point[] = {0, 0, 0};

    boost::array<double, N_pixel * N_pixel> data;
    data.fill(0.0);

    // const double cutoff(20000);
    const double cutoff(2000);
    double Itot(0.0);

    // const unsigned int N_point(17200);
    // const unsigned int N_point(1000);
    const unsigned int N_point(2620);

    if (argc == 1)
    {
        double points[N_point][3];
        double intensity[N_point];
        generate_random_points(points, N_point, L);
        emission(points, intensity, N_point);

        for (unsigned int i(0); i < N_point; ++i)
        {
            std::cout << "[" << i << "/"
                << N_point << "] = " << point_as_str(points[i]) << std::endl;

            overlay_psf(data.data(), N_pixel, pixel_length,
                points[i], intensity[i], focal_point, k, N_A, cutoff);
            Itot += intensity[i];
        }
    }
    else
    {
        double shift[3] = {15, 15, 1.5};
        double scale = 1000.0;

        const unsigned int frames(argc - 1);

        #ifdef _OPENMP
        std::vector<boost::array<double, N_pixel * N_pixel> > tmp(frames);

        #pragma omp parallel for reduction(+:Itot)
        #endif
        for (unsigned int cnt = 0; cnt < frames; ++cnt)
        {
            // std::string const filename(
            //     (boost::format("sample_data3/test%03d.csv")
            //         % (cnt * (101 / frames))).str());
            std::string const filename(argv[cnt + 1]);

            double points[N_point][3];
            double intensity[N_point];
            read_input(filename.c_str(), points, N_point, shift, scale);
            emission(points, intensity, N_point);

            #ifdef _OPENMP
            tmp[cnt].fill(0.0);
            #endif

            for (unsigned int i(0); i < N_point; ++i)
            {
                std::cout << "[" << cnt << ":" << i << "/"
                    << N_point << "] = " << point_as_str(points[i]) << std::endl;

                const double I(intensity[i] / frames);
                #ifdef _OPENMP
                overlay_psf(tmp[cnt].data(), N_pixel, pixel_length,
                    points[i], I, focal_point, k, N_A, cutoff);
                #else
                overlay_psf(data.data(), N_pixel, pixel_length,
                    points[i], I, focal_point, k, N_A, cutoff);
                #endif
                Itot += I;
            }
        }

        #ifdef _OPENMP
        for (unsigned int cnt = 0; cnt < frames; ++cnt)
        {
            for (unsigned int i = 0; i < N_pixel * N_pixel; ++i)
            {
                data[i] += tmp[cnt][i];
            }
        }
        #endif
    }

    const double I(boost::accumulate(data, 0.0));
    std::cout << "The total intensity of an output image is " << I << std::endl;
    std::cout << "The total intensity of PSFs is " << Itot << std::endl;
    std::cout << "The relative error is " << fabs(Itot - I) / Itot << std::endl;

    detection(data.data(), data.data(), N_pixel * N_pixel);
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
    save_data("result.txt", data.data(), N_pixel * N_pixel);
}

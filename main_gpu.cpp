#include "microscope_gpu.h"

extern "C" void born_wolf_psf_tbl_gpu(
    double data[], unsigned int const N_pixel, double const pixel_length,
    double points[][3], unsigned int const N_point,
    double const I, double c[3], double const k, double const N_A,
    double const cutoff);

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

#include <boost/array.hpp>

std::string point_as_str(double p[3])
{
    std::stringstream sout;
    sout << std::showpos;
    sout << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return sout.str();
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
            // t = val;
            ss.str("");
            ss.clear(std::stringstream::goodbit);
        }
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

        points[i][0] = (y - shift[0]) * scale;
        points[i][1] = (z - shift[1]) * scale;
        points[i][2] = fabs(x - shift[2]) * scale;
        // intensity[i] = emission(points[i][2]);
        // std::cout << point_as_str(points[i]) << std::endl;
        i++;
    };
}

int main(int argc, char** argv)
{
    const unsigned int N_point(1000);
    double points[N_point][3];
    double shift[3] = {15e-6, 15e-6, 5.06228e-07};
    double scale = 1000.0 * 1e+6;

    const double lambda = 508.0; // nm
    const double N_A = 1.4;
    const double k = 2 * M_PI / lambda;
    double focal_point[] = {0, 0, 0};
    const double cutoff = 2000.0;
    const unsigned int N_pixel = 600;
    const double objective = 100.0;
    const double pixel_length = 6500 / objective;

    boost::array<double, N_pixel * N_pixel> data;
    data.fill(0.0);

    for (unsigned int idx = 1; idx < argc; ++idx)
    {
        const std::string filename(argv[idx]);
        read_input(filename.c_str(), points, N_point, shift, scale);

        const double intensity = 1.0;

        born_wolf_psf_tbl_gpu(
            data.data(), N_pixel, pixel_length, points, N_point,
            intensity, focal_point, k, N_A, cutoff);
    }

    std::cout << std::scientific << std::setprecision(16);
    for (unsigned int i = 0; i < N_pixel * N_pixel; ++i)
    {
        std::cout << data[i] << std::endl;
    }

    return 0;
}

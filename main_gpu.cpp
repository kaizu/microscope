#include "microscope_gpu.h"

extern "C" void born_wolf_psf_tbl_gpu(
    double data[], unsigned int const N_pixel, double const pixel_length,
    double points[][3], unsigned int const N_point,
    double intensity[], double c[3], double const k, double const N_A,
    double const cutoff);

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

#include <boost/array.hpp>

#include "microscope.hpp"

using namespace microscope;

std::string point_as_str(double p[3])
{
    std::stringstream sout;
    sout << std::showpos;
    sout << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
    return sout.str();
}

void read_input(
    char const filename[], double points[][3], double intensity[], unsigned int const data_size)
{
    std::ifstream fin(filename);
    std::string buf;

    // std::getline(fin, buf); // Ignore a header line

    if (fin)
    {
        unsigned int i(0);
        while (std::getline(fin, buf) && i < data_size)
        {
            std::string token;
            std::istringstream stream(buf);
            std::stringstream ss;
            double val;
            unsigned int sid;

            double x, y, z, I;

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
                ss << token;
                ss >> val;
                I = val;
                ss.str("");
                ss.clear(std::stringstream::goodbit);
            }

            points[i][0] = x;
            points[i][1] = y;
            points[i][2] = z;
            intensity[i] = I;
            // intensity[i] = emission(points[i][2]);
            // std::cout << point_as_str(points[i]) << std::endl;
            i++;
        };
    }
}

std::string modify_filename(
    const std::string& filename, const std::string& prefix, const std::string& suffix)
{
    std::string dirname, basename;

    const size_t lastsep = filename.find_last_of("/");
    if (lastsep != std::string::npos)
    {
        dirname = filename.substr(0, lastsep);
        basename = filename.substr(lastsep + 1, filename.size());
    }
    else
    {
        dirname = ".";
        basename = filename;
    }

    std::string root, ext;

    const size_t lastdot = basename.find_last_of(".");
    if (lastdot != std::string::npos)
    {
        if (lastdot == 0)
        {
            root = basename;
            ext = "";
        }
        else
        {
            root = basename.substr(0, lastdot);
            ext = basename.substr(lastdot, basename.size());
        }
    }
    else
    {
        root = basename;
        ext = "";
    }

    const std::string outname = dirname + "/" + prefix + root + suffix;
    // const std::string outname = root + "." + newext;
    std::cout << outname << std::endl;
    return outname;
}

bool save_data(char const filename[], double data[], unsigned int const data_size)
{
    std::ofstream fout;
    fout.open(filename);
    if (fout.fail())
    {
        return false;
    }
    fout.setf(std::ios::scientific);
    fout.precision(16);
    for (unsigned int i(0); i < data_size; ++i)
    {
        fout << data[i] << std::endl;
    }
    fout.close();
    return true;
}

int main(int argc, char** argv)
{
    const unsigned int N_point(1000);
    double points[N_point][3];

    // const double wave_length = 508.0;
    // const double N_A = 1.4;
    // const double k = 2 * M_PI / wave_length;
    // const double cutoff = 2000.0;
    // const unsigned int N_pixel = 600;
    // const double objective = 100.0;
    // const double pixel_length = 6500 / objective;
    // const double angle = 75.0;
    // const double flux_density = 20.0;
    // const std::pair<double, double> field
    //     = tirf_evanescent_field(angle, wave_length, flux_density);
    // const double ATsq = field.first;
    // const double d = field.second;
    // // const double epsilon = 0.0145;
    // const double epsilon = 0.00435;
    // const double absorption_cross_section = 3.19e-18;
    // const double exposure_time = 30e-3;

    // const double alpha = N_A * k;
    // double focal_point[] = {0, 0, 0};
    // // double focal_point[] = {0, 0, -240.0};

    const double wave_length = 512.0;
    const double N_A = 1.49;
    const double k = 2 * M_PI / wave_length;
    const double cutoff = 2000.0;
    const unsigned int N_pixel = 512;
    const double objective = 60.0 * 1.5 * 4.0;
    // const double pixel_length = 16000 / objective;
    const double pixel_length = 44.0;
    const double angle = 75.0;
    const double flux_density = 20.0;
    const std::pair<double, double> field
        = tirf_evanescent_field(angle, wave_length, flux_density);
    const double ATsq = field.first;
    const double d = field.second;
    const double epsilon = 0.5 / 300;
    // const double epsilon = 1.5 * 2.0 / 300;
    // const double epsilon = 0.3;
    // const double epsilon = 0.0145;
    // const double epsilon = 0.00435;
    const double absorption_cross_section = 3.19e-18;
    const double exposure_time = 30e-3;

    const double alpha = N_A * k;
    double focal_point[] = {0, 0, 0};
    // double focal_point[] = {0, 0, -5 * pixel_length};
    // double focal_point[] = {0, 0, -240.0};

    boost::array<double, N_pixel * N_pixel> data;
    boost::array<double, N_point> intensity;

    for (unsigned int idx = 1; idx < argc; ++idx)
    {
        const std::string filename(argv[idx]);
        read_input(filename.c_str(), points, intensity.data(), N_point);

        // intensity.fill(1.0);
        emission(
            points, intensity.data(), N_point,
            ATsq, d, epsilon, absorption_cross_section, exposure_time);

        data.fill(0.0);
        born_wolf_psf_tbl_gpu(
            data.data(), N_pixel, pixel_length, points, N_point,
            intensity.data(), focal_point, k, N_A, cutoff);

        emccd_detection(data.data(), data.data(), data.size());

        save_data(modify_filename(filename, "setting4_", ".dat").c_str(), data.data(), data.size());
        // save_data(modify_filename(filename, "setting2_", ".dat").c_str(), data.data(), data.size());
    }

    return 0;
}

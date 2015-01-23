#include <iostream>
#include <fstream>
#include <sstream>
#include <gsl/gsl_rng.h>

#include <boost/array.hpp>
#include <boost/range/numeric.hpp>

// #include "point_spreading_functions.hpp"
#include "born_wolf_psf.hpp"
using namespace microscope;


int main()
{
    const unsigned int N(2000 + 1);
    const unsigned int M(200 + 1);
    const double rmax(675);
    const double zmax(60);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(16);

    const double dr(rmax / N);
    const double dz(zmax / M);

    std::ofstream fout;
    fout.open("born_wolf_psf_table.hpp");
    fout.setf(std::ios::scientific);
    fout.precision(16);

    fout << "#ifndef __BORN_WOLF_PSF_TABLE_HPP" << std::endl;
    fout << "#define __BORN_WOLF_PSF_TABLE_HPP" << std::endl;
    fout << std::endl;
    fout << "namespace born_wolf_psf_table" << std::endl;
    fout << "{" << std::endl;
    fout << std::endl;

    fout << "struct Table" << std::endl;
    fout << "{" << std::endl;
    fout << "    const unsigned int N;" << std::endl;
    fout << "    const unsigned int M;" << std::endl;
    fout << "    const double rmax;" << std::endl;
    fout << "    const double zmax;" << std::endl;
    fout << "    const double* const y;" << std::endl;
    fout << "};" << std::endl;
    fout << std::endl;

    fout << "static const double born_wolf_psf_table_f"
        << "[" << M  * N << "] = {" << std::endl;
    for (unsigned int m(0); m < M; ++m)
    {
        const double z(m * dz);

        for (unsigned int n(0); n < N; ++n)
        {
            const double r(n * dr);
            const double val(born_wolf_psf(r, z, 1.0, 1.0));

            fout << "        " << val
                << (n != N - 1 || m != M - 1 ? "," : "") << std::endl;
        }

        std::cout << m << std::endl;
    }
    fout << "};" << std::endl;
    fout << std::endl;

    // fout << "static const double born_wolf_psf_table_f"
    //     << "[" << M - 1 << " + 1][" << N - 1 << " + 1] = {" << std::endl;
    // for (unsigned int m(0); m < M; ++m)
    // {
    //     const double z(m * dz);

    //     fout << "    {" << std::endl;

    //     double tot(0.0);
    //     for (unsigned int n(0); n < N; ++n)
    //     {
    //         const double r(n * dr);
    //         const double val(born_wolf_psf(r, z, 1.0, 1.0));
    //         tot += val * r * (rmax / N);

    //         fout << "        " << val << (n != N - 1 ? "," : "") << std::endl;
    //     }

    //     std::cout << m << "," << tot << ","
    //         << int_psf_cylinder(0, rmax, z, 1.0, 1.0) << std::endl;

    //     fout << "    }" << (m != M - 1 ? "," : "") << " // " << m << std::endl;
    // }
    // fout << "};" << std::endl;
    // fout << std::endl;

    fout << "static const Table born_wolf_psf_table = {"
        << N << ", " << M << ", " << rmax << ", " << zmax
        << ", born_wolf_psf_table_f};" << std::endl;

    fout << "} // born_wolf_psf_table" << std::endl;
    fout << std::endl;
    fout << "#endif /* __BORN_WOLF_PSF_TABLE_HPP */" << std::endl;

    fout.close();
}

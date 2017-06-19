#include "point_spreading_functions.hpp"

int main(void)
{
    const double wave_length = 508.0;
    const double N_A = 1.4;
    const double k = 2 * M_PI / wave_length;

    const unsigned int N = born_wolf_psf_table::born_wolf_psf_table.N;
    const double rmax = born_wolf_psf_table::born_wolf_psf_table.rmax;
    const unsigned int M = born_wolf_psf_table::born_wolf_psf_table.M;
    const double zmax = born_wolf_psf_table::born_wolf_psf_table.zmax;

    for (unsigned int i = 0; i < N; ++i)
    {
        const double r = 0.0 + rmax * (double)(i) / (double)(N);
        std::cout << r;
        // const double z = 0.0 + zmax * (double)(2) * 100.0 / (double)(M);
        // std::cout << " " << microscope::born_wolf_psf(r, z, k, N_A);
        // std::cout << " " << microscope::born_wolf_psf_tbl(r, z, k, N_A);
        for (unsigned int j = 0; j < 5; ++j)
        {
            const double z = 0.0 + zmax * (double)(j) * 100.0 / (double)(M);
            const double I = microscope::born_wolf_psf_tbl(r, z, k, N_A);
            // const double I = microscope::born_wolf_psf(r, z, k, N_A);
            std::cout << " " << I;
        }
        std::cout << std::endl;
    }
    return 0;
}

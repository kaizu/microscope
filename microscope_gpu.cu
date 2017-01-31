#include "microscope_gpu.h"

extern "C" void born_wolf_psf_tbl_gpu(
    double data[], unsigned int const N_pixel, double const pixel_length,
    double points[][3], unsigned int const N_point,
    double const I, double c[3], double const k, double const N_A,
    double const cutoff);

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include "born_wolf_psf_table.hpp"

#define PIXELS 63
#define BLOCKS 63 * 63
#define THREADS 128
#define N_DIV 11

typedef struct
{
    double px;
    double py;
    double xmin;
    double ymin;
    unsigned int N;
    unsigned int M;
    double imin;
    double jmin;
    unsigned int N_pixel;
    double pixel_length;
    double delta_length;
    double alpha_factor;
    double psi_factor;
    double C;
} parameter_type;

template <unsigned int BLOCK_SIZE>
__global__ void born_wolf_psf_tbl_gpu_kernel(
    double* d_y, double* d_data, const parameter_type params)
{
    extern __shared__ double s_x[];
    const unsigned int tid = threadIdx.x;
    // const unsigned int bid = blockIdx.y * gridDim.x + blockIdx.x;
	// const unsigned int idx = bid * blockDim.x + tid;

    if (tid < N_DIV * N_DIV)
    {
        const double xmin = params.xmin + params.pixel_length * blockIdx.x;
        const double ymin = params.ymin + params.pixel_length * blockIdx.y;
        // const double xmax = xmin + params.pixel_length;
        // const double ymax = ymin + params.pixel_length;

        const unsigned int xidx = tid % N_DIV;
        const unsigned int yidx = tid / N_DIV;

        const double x = xmin + params.delta_length * xidx;
        const double y = ymin + params.delta_length * yidx;
        const double dx = x - params.px;
        const double dy = y - params.py;
        const double r = sqrt(dx * dx + dy * dy);

        const double n1 = r * params.alpha_factor;
        const double n2 = floor(n1);
        const double m1 = 0.0 * params.psi_factor;
        const double m2 = floor(m1);
        const unsigned int n = static_cast<unsigned int>(n2);
        const unsigned int m = static_cast<unsigned int>(m2);

        if (n < params.N && m < params.M)
        {
            const double pr = n1 - n2;
            const double pz = m1 - m2;
            const double v00 = d_y[m * params.N + n];
            const double v10 = d_y[m * params.N + n + 1];
            const double v01 = d_y[(m + 1) * params.N + n];
            const double v11 = d_y[(m + 1) * params.N + n + 1];
            const double I = (1.0 - pr) * ((1.0 - pz) * v00 + pz * v01) + pr * ((1.0 - pz) * v10 + pz * v11);

            const double factor1 = ((xidx == 0 || xidx == N_DIV - 1) ? 1.0 / 3.0 : (xidx % 2 == 0 ? 2.0 / 3.0 : 4.0 / 3.0));
            const double factor2 = ((yidx == 0 || yidx == N_DIV - 1) ? 1.0 / 3.0 : (yidx % 2 == 0 ? 2.0 / 3.0 : 4.0 / 3.0));

            s_x[tid] = factor1 * factor2 * params.delta_length * params.delta_length * I;
        }
    }

    __syncthreads();

    // for (unsigned int i = 1; i < blockDim.x; i *= 2)
    // {
    //     if (tid % (2 * i) == 0)
    //     {
    //         s_x[tid] += s_x[tid + i];
    //     }
    //     __syncthreads();
    // }

    // for (unsigned int i = blockDim.x / 2; i > 0; i >>= 1)
    // {
    //     if (tid < i)
    //     {
    //         s_x[tid] += s_x[tid + i];
    //     }
    //     __syncthreads();
    // }

    if (BLOCK_SIZE >= 512)
    {
        if (tid < 256)
        {
            s_x[tid] += s_x[tid + 256];
        }
        __syncthreads();
    }

    if (BLOCK_SIZE >= 256)
    {
        if (tid < 128)
        {
            s_x[tid] += s_x[tid + 128];
        }
        __syncthreads();
    }

    if (BLOCK_SIZE >= 128)
    {
        if (tid < 64)
        {
            s_x[tid] += s_x[tid + 64];
        }
        __syncthreads();
    }

    if (tid < 32)
    {
        if (BLOCK_SIZE >= 64) s_x[tid] += s_x[tid + 32];
        if (BLOCK_SIZE >= 32) if (tid < 16) s_x[tid] += s_x[tid + 16];
        if (BLOCK_SIZE >= 16) if (tid < 8) s_x[tid] += s_x[tid + 8];
        if (BLOCK_SIZE >= 8) if (tid < 4) s_x[tid] += s_x[tid + 4];
        if (BLOCK_SIZE >= 4) if (tid < 2) s_x[tid] += s_x[tid + 2];
        if (BLOCK_SIZE >= 2) if (tid < 1) s_x[tid] += s_x[tid + 1];
    }

    if (tid == 0)
    {
        const int i1 = static_cast<int>(static_cast<double>(blockIdx.x) + params.imin);
        const int j1 = static_cast<int>(static_cast<double>(blockIdx.y) + params.jmin);
        if (0 <= i1 && i1 < params.N_pixel && 0 <= j1 && j1 < params.N_pixel)
        {
            d_data[static_cast<unsigned int>(j1 * params.N_pixel + i1)] += s_x[0] * params.C;
        }
    }
}

/*
 * const unsigned int N_pixel = 600;
 * const double objective = 100.0;
 * const double pixel_length = 6500.0 / objective;
 * const double k = 2.0 * M_PI / 508.0;
 * const double N_A = 1.4;
 * const unsigned int N_point = 1000;
 */

void born_wolf_psf_tbl_gpu(
    double data[], unsigned int const N_pixel, double const pixel_length,
    double points[][3], unsigned int const N_point,
    double const I, double c[3], double const k, double const N_A,
    double const cutoff)
{
    const unsigned int N = born_wolf_psf_table::born_wolf_psf_table.N;
    const unsigned int M = born_wolf_psf_table::born_wolf_psf_table.M;
    const double rmax = born_wolf_psf_table::born_wolf_psf_table.rmax;
    const double zmax = born_wolf_psf_table::born_wolf_psf_table.zmax;

    const double alpha = N_A * k;
    const double C = alpha * alpha / M_PI;
    const double alpha_factor = alpha * N / rmax;
    const double psi_factor = 0.5 * alpha * N_A * M / zmax;
    const double offset = N_pixel * pixel_length * -0.5;
    const double delta_length = pixel_length / (N_DIV - 1);

    const unsigned int tbl_size = N * M; // 802401;
    thrust::host_vector<double> y(tbl_size);
    memcpy((void *)thrust::raw_pointer_cast(y.data()),
           (const void *)born_wolf_psf_table::born_wolf_psf_table.y,
           tbl_size * sizeof(double));
    thrust::device_vector<double> d_y = y;

    thrust::device_vector<double> d_x(N_pixel * N_pixel, 0.0);
    // thrust::host_vector<double> x(N_pixel * N_pixel, 0.0);

    for (unsigned int cnt = 0; cnt < N_point; ++cnt)
    {
        const double x = points[cnt][0] - c[0];
        const double y = points[cnt][1] - c[1];

        const double i0 = floor((x - offset) / pixel_length);
        const double j0 = floor((y - offset) / pixel_length);
        const double imin = i0 - floor(PIXELS * 0.5);
        const double jmin = j0 - floor(PIXELS * 0.5);
        const double xmin = imin * pixel_length + offset;
        const double ymin = jmin * pixel_length + offset;

        parameter_type params = {
            x, y, xmin, ymin,
            N, M, imin, jmin, N_pixel,
            pixel_length, delta_length,
            alpha_factor, psi_factor, C};

        born_wolf_psf_tbl_gpu_kernel<THREADS><<<dim3(PIXELS, PIXELS, 1), THREADS, THREADS * sizeof(double)>>>(thrust::raw_pointer_cast(d_y.data()), thrust::raw_pointer_cast(d_x.data()), params);
    }

    // x = d_x;
    cudaMemcpy(data, thrust::raw_pointer_cast(d_x.data()), (N_pixel * N_pixel) * sizeof(double), cudaMemcpyDeviceToHost);
}

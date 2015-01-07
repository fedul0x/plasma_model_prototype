#include <math.h>
__device__ __host__ double dif(float x, float y, float z, float step)
{
    return (double)((x - y - y + z) / step / step);
}

//TODO change step and matrix_size to vector
__global__ void potential_establish(float *prev_phi, float *next_phi, float *flags, float *step, int *matrix_size )
{
    const int x_dim = matrix_size[0];
    const int y_dim = matrix_size[1];
    const int z_dim = matrix_size[2];

    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (x_dim - 2) * (y_dim - 2) * (z_dim - 2))
        return;
    const int i = (int) (index / (z_dim - 2) / (y_dim - 2)) + 1;
    const int j = ((int) (index / (z_dim - 2))) % (y_dim - 2) + 1;
    const int k = index % (z_dim - 2) + 1;

    const float x_step = step[0];
    const float y_step = step[1];
    const float z_step = step[2];

    const int ijk = i * z_dim * y_dim + j * z_dim + k;
    const int im1jk = (i - 1) * z_dim * y_dim + j * z_dim + k;
    const int ijm1k = i * z_dim * y_dim + (j - 1) * z_dim + k;
    const int ijkm1 = i * z_dim * y_dim + j * z_dim + k - 1;
    const int ip1jk = (i + 1) * z_dim * y_dim + j * z_dim + k;
    const int ijp1k = i * z_dim * y_dim + (j + 1) * z_dim + k;
    const int ijkp1 = i * z_dim * y_dim + j * z_dim + k + 1;

    const int indexes[6] = {im1jk, ijm1k, ijkm1, ip1jk, ijp1k, ijkp1};

    float tmp = prev_phi[ijk];
    float ro = 1E-17;
    float sum;
    int cc = 0;
    for (int count = 1; count <= 200; count++)
    {
        prev_phi[ijk] = tmp;
        flags[ijk] = count;

        // do
        // {
        //     sum = 0;
        //     cc = 0;
        //     for (int b = 0; b < 6; b++)
        //     {
        //         if (indexes[b] != 0)
        //         {
        //             cc += 1;
        //             sum += flags[indexes[b]];
        //         }

        //     }

        // }  while (sum/cc == count-1);
        // tmp = 1.75E-10 * (
        //           (prev_phi[ip1jk] - 2*prev_phi[ijk] + prev_phi[im1jk]) / x_step / x_step + 
        //           (prev_phi[ijp1k] - 2*prev_phi[ijk] + prev_phi[ijm1k]) / y_step / y_step +
        //           (prev_phi[ijkp1] - 2*prev_phi[ijk] + prev_phi[ijkm1]) / z_step / z_step + 
        //           // dif(prev_phi[ip1jk], prev_phi[ijk], prev_phi[im1jk], x_step) +
        //           // dif(prev_phi[ijp1k], prev_phi[ijk], prev_phi[ijm1k], y_step) +
        //           // dif(prev_phi[ijkp1], prev_phi[ijk], prev_phi[ijkm1], z_step) +
        //            4 * M_PI * ro) + prev_phi[ijk];
        tmp = 1.75E-10 * (
                  (prev_phi[ip1jk] - 2*prev_phi[ijk] + prev_phi[ijk]) / x_step / x_step + 
                  (prev_phi[ijk] - 2*prev_phi[ijk] + prev_phi[ijk]) / y_step / y_step +
                  (prev_phi[ijk] - 2*prev_phi[ijk] + prev_phi[ijk]) / z_step / z_step + 
                   4 * M_PI * ro) + prev_phi[ijk];
        //            __syncthreads();
        // next_phi[ijk] = tmp;
        next_phi[ijk] = tmp;
        // prev_phi[ijk] = tmp;
        // index += 512;
    }
}
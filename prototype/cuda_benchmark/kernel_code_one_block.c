#include <math.h>
__device__ __host__ double dif(float x, float y, float z, float step)
{
    return (double)((x - y - y + z) / step / step);
}

//TODO change step and matrix_size to vector
__global__ void potential_establish(float *prev_phi, float *next_phi, float *step, int *matrix_size )
{
    const int x_dim = matrix_size[0];
    const int y_dim = matrix_size[1];
    const int z_dim = matrix_size[2];

    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= (x_dim - 2) * (y_dim - 2) * (z_dim - 2))
        return;
    int i = (int) (index / (z_dim - 2) / (y_dim - 2)) + 1;
    int j = ((int) (index / (z_dim - 2))) % (y_dim - 2) + 1;
    int k = index % (z_dim - 2) + 1;

    const float x_step = step[0];
    const float y_step = step[1];
    const float z_step = step[2];
    // __shared__ int indexes[7];
    

    __shared__ int ijk;
    __shared__ int im1jk;
    __shared__ int ijm1k;
    __shared__ int ijkm1;
    //__shared__ int ijkm1;
    __shared__ int ip1jk;
    __shared__ int ijp1k;
    __shared__ int ijkp1;
    //__shared__ int ijkp1;
    float tmp = prev_phi[ijk];
    float ro = 1E-17;
    float precision = 0;
    const int limit = 0;//(x_dim - 2);
    //for (int l = 0; l <= 5; l++)
    {
        index = blockIdx.x * blockDim.x + threadIdx.x;
        for (int l = 0; l<2; l++)
         // while (index < limit)
        {
            i = (int) (index / (z_dim - 2) / (y_dim - 2)) + 1;
            j = ((int) (index / (z_dim - 2))) % (y_dim - 2) + 1;
            k = index % (z_dim - 2) + 1;

            ijk = i * z_dim * y_dim + j * z_dim + k;
            im1jk = (i - 1) * z_dim * y_dim + j * z_dim + k;
            ijm1k = i * z_dim * y_dim + (j - 1) * z_dim + k;
            ijkm1 = i * z_dim * y_dim + j * z_dim + k - 1;
            ip1jk = (i + 1) * z_dim * y_dim + j * z_dim + k;
            ijp1k = i * z_dim * y_dim + (j + 1) * z_dim + k;
            ijkp1 = i * z_dim * y_dim + j * z_dim + k + 1;

            // tmp = prev_phi[ijk];
            // prev_phi[ijk] = tmp;
            precision = (prev_phi[ip1jk] - prev_phi[ijk] - prev_phi[ijk] + prev_phi[im1jk]) / x_step / x_step +
                        (prev_phi[ijp1k] - prev_phi[ijk] - prev_phi[ijk] + prev_phi[ijm1k]) / y_step / y_step +
                        (prev_phi[ijkp1] - prev_phi[ijk] - prev_phi[ijk] + prev_phi[ijkm1]) / z_step / z_step +
                        4 * M_PI * ro;
            tmp = (float)(1.75E-10 * (precision) + prev_phi[ijk]);
            // tmp = 1;
            __syncthreads();
            next_phi[ijk] = tmp;
            prev_phi[ijk] = tmp;
            // index += 512;
        }
    }
}
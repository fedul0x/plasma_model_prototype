#include <math.h>
__device__ __host__ float dif(float x, float y, float z, float step){
    return (x - 2 * y + z) / step / step;
}

//TODO change step and matrix_size to vector
__global__ void potential_establish(float* prev_phi, float* next_phi, float* step, int* matrix_size )
{
    const int x_dim = matrix_size[0];
    const int y_dim = matrix_size[1];
    const int z_dim = matrix_size[2];

    const int index = blockIdx.x*blockDim.x + threadIdx.x;
    if (index >=(x_dim-2)*(y_dim-2)*(z_dim-2)) 
        return;
    const int i = (int) (index/(z_dim-2)/(y_dim-2))+1;
    const int j = ((int) (index/(z_dim-2)))%(y_dim-2)+1;
    const int k = index%(z_dim-2)+1;
    // const int i = (int) (threadIdx.x/(z_dim-2)/(y_dim-2))+1;
    // const int j = ((int) (threadIdx.x/(z_dim-2)))%(y_dim-2)+1;
    // const int k = threadIdx.x%(z_dim-2)+1;

    const float x_step = step[0];
    const float y_step = step[1];
    const float z_step = step[2];

    const int ijk = i*z_dim*y_dim + j*z_dim + k;
    const int im1jk = (i-1)*z_dim*y_dim + j*z_dim + k;
    const int ijm1k = i*z_dim*y_dim + (j-1)*z_dim + k;
    const int ijkm1 = i*z_dim*y_dim + j*z_dim + k - 1;
    //const int ijkm1 = ijk - 1;
    const int ip1jk = (i+1)*z_dim*y_dim + j*z_dim + k;
    const int ijp1k = i*z_dim*y_dim + (j+1)*z_dim + k;
    const int ijkp1 = i*z_dim*y_dim + j*z_dim + k + 1;
    //const int ijkp1 = ijk + 1;
    float tmp = prev_phi[ijk];
    float ro = 1E-17;
    for(int l=0; l <= 500; l++) 
    {
        prev_phi[ijk] = tmp;
        //TODO rewrite PI and deltaT_
        //next_phi[i][j][k] = 1.75E-10 * (dif(prev_phi[i + 1][j][k], prev_phi[i][j][k], prev_phi[i - 1][j][k], x_step) + dif(prev_phi[i][j + 1][k], prev_phi[i][j][k], prev_phi[i][j - 1][k], y_step) + dif(prev_phi[i][j][k + 1], prev_phi[i][j][k], prev_phi[i][j][k - 1], z_step) + 4 * 3.1415926 * ro[i][j][k]) + prev_phi[i][j][k];
        tmp = 1.75E-10 * (
        dif(prev_phi[ip1jk], prev_phi[ijk], prev_phi[im1jk], x_step) +
        dif(prev_phi[ijp1k], prev_phi[ijk], prev_phi[ijm1k], y_step) +
        dif(prev_phi[ijkp1], prev_phi[ijk], prev_phi[ijkm1], z_step) +
        4 * M_PI * ro) + prev_phi[ijk];
        __syncthreads();
        //TODO change this append to another where not depend to next 
        next_phi[ijk] = tmp;
        //next_phi[ijk] = next_phi[ijk] + 1;
    }
}
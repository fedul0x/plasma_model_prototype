#include <math.h>


__device__ volatile unsigned int count;                 // Флаг-счетчик отработавших блоков. Под него выделится
__device__ volatile unsigned int term;                 // Флаг-счетчик отработавших блоков. Под него выделится
//4 байта в глобальной памяти устройства.

/* Фунция начальной инициализации флага-счетчика: */
__device__ void InitSyncWholeDevice(const int index)
{
    if (index == 0)                            // Первый поток в grid`е (индекс 0) запишет нулевое
        count = 0;                             //начальное значение в счетчик блоков.
    term = 5000;
    // term--;
    if (threadIdx.x == 0)                      // Первый поток каждого block`а будет ждать, пока флаг-
        while (count != 0);                    //счетчик действительно станет нулем.
    // while ((count != 0) && (term >= 0))                    //счетчик действительно станет нулем.
    //     {
    //         atomicDec((unsigned int*)&term, gridDim.x-1);
    //     }

    // Заставляем остальные потоки каждого block`а ждать, пока первые не выйдут из цикла:
    __syncthreads();
    // Все, флаг-аккумулятор записан. Все потоки на device более-менее идут вровень.
}

/* Фунция синхронизации потоков на device: */
__device__ void SyncWholeDevice()
{
    // Переменная под значение счетчика до инкремента:
    unsigned int oldc;
    // Каждый поток пождет, пока записанное им в gmem и smem, станет видно всему grid`у:
    term = 1;
    __threadfence();
    term = 2;

    // Первые потоки каждого block`а атомарным образом инкрементируют (каждый по разу)
    //флаг-аккумулятор:
    if (threadIdx.x == 0)
    {
        term = 3;
        // В oldc кладется значение count до "+1":
        // oldc = atomicInc(&count, gridDim.x-1);
        // oldc = atomicInc((unsigned int *)&count, gridDim.x - 1);
        oldc = atomicInc((unsigned int *)&count, gridDim.x-1);
        term = 4;
        // Пусть поток подождет, пока его инкремент "дойдет" до ячейки в gmem:
        __threadfence();
        term = 5;

        // Если это последний блок (остальные уже инкрементировали count и ждут за счет цикла ниже),
        //то и незачем ему считывать count, так как предварительно убедились, что его инкремент
        //записан в gmem. Если мы в блоке, который еще не "отработал", то его первый поток будет
        //зациклен, пока все остальные блоки не "отчитаются" о завершении счета.
        // long term = 16000;
        if (oldc != (gridDim.x - 1))
            while (count != 0);
            // while ((count != 0) && (term-- != 0));
            // while (count != gridDim.x);
    }
    term = 6;
    // return;

    // Заставляем потоки в каждом блоке ждать, пока первые не выйдут из цикла:
    __syncthreads();
}


__device__ __host__ float dif(float x, float y, float z, float step)
{
    return (x - 2 * y + z) / step / step;
}

//TODO change step and matrix_size to vector
__global__ void potential_establish(float *prev_phi, float *next_phi, float *step, int *matrix_size )
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
    // const int i = (int) (threadIdx.x/(z_dim-2)/(y_dim-2))+1;
    // const int j = ((int) (threadIdx.x/(z_dim-2)))%(y_dim-2)+1;
    // const int k = threadIdx.x%(z_dim-2)+1;

    const float x_step = step[0];
    const float y_step = step[1];
    const float z_step = step[2];

    const int ijk = i * z_dim * y_dim + j * z_dim + k;
    const int im1jk = (i - 1) * z_dim * y_dim + j * z_dim + k;
    const int ijm1k = i * z_dim * y_dim + (j - 1) * z_dim + k;
    const int ijkm1 = i * z_dim * y_dim + j * z_dim + k - 1;
    //const int ijkm1 = ijk - 1;
    const int ip1jk = (i + 1) * z_dim * y_dim + j * z_dim + k;
    const int ijp1k = i * z_dim * y_dim + (j + 1) * z_dim + k;
    const int ijkp1 = i * z_dim * y_dim + j * z_dim + k + 1;
    //const int ijkp1 = ijk + 1;
    float tmp = prev_phi[ijk];
    float ro = 1E-17;
    for (int l = 0; l <= 5; l++)
    {
        InitSyncWholeDevice(threadIdx.x + blockIdx.x * blockDim.x);
        // return;
        prev_phi[ijk] = tmp;
        //TODO rewrite PI and deltaT_
        //next_phi[i][j][k] = 1.75E-10 * (dif(prev_phi[i + 1][j][k], prev_phi[i][j][k], prev_phi[i - 1][j][k], x_step) + dif(prev_phi[i][j + 1][k], prev_phi[i][j][k], prev_phi[i][j - 1][k], y_step) + dif(prev_phi[i][j][k + 1], prev_phi[i][j][k], prev_phi[i][j][k - 1], z_step) + 4 * 3.1415926 * ro[i][j][k]) + prev_phi[i][j][k];
        tmp = 1.75E-10 *
              (dif(prev_phi[ip1jk], prev_phi[ijk], prev_phi[im1jk], x_step) +
               dif(prev_phi[ijp1k], prev_phi[ijk], prev_phi[ijm1k], y_step) +
               dif(prev_phi[ijkp1], prev_phi[ijk], prev_phi[ijkm1], z_step) +
               4 * M_PI * ro) + prev_phi[ijk];
        //TODO change this append to another where not depend to next
        next_phi[ijk] = tmp;
        prev_phi[ijk] = tmp;
        //next_phi[ijk] = next_phi[ijk] + 1;
        SyncWholeDevice();
        // next_phi[0] = count;
        // if (l == 0) 
        return;
    }
}
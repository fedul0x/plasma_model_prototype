#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define positionElectron pe


struct Point
{
    double x;
    double y;
    double z;
};

/* Перевод скорости к безразмерному виду */
long double convertSpeedToDimensionless(long double speed, long double Mve)
{
    return speed / Mve;
}

/* Вычисление расстояния (модуля вектора) */
long double calcDistance(struct Point p)
{
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

long double calcDistance(long double x, long double y, long double z)
{
    return sqrt(x * x + y * y + z * z);
}

double randFrom0To1()
{
    return (double)rand() / (double)RAND_MAX;
}

int main(int argc, char const *argv[])
{
    unsigned int r = 2000000;
    srand(r);

    /* Число частиц в ячейке */
    long double numbersElectron = 1E+10;
    long double numbersCarbon = 1E+16;
    long double numbersHelium = 1E+10;
    /* Массы частиц */
    long double massElectron = 9.11E-31;
    long double massCarbon = 12.011 * 1.67E-27 - 9.11E-31;
    long double massHelium = 4.002 * 1.67E-27 - 9.11E-31;
    /* Заряды частиц */
    long double chargeElectron = -1.602176565E-19;
    long double chargeCarbon = 1.602176565E-19;
    long double chargeHelium = 1.602176565E-19;
    /* Радиусы частиц */
    long double radiusElectron = 2.81794E-15;
    long double radiusCarbon = 31E-12;
    long double radiusHelium = 91E-12;
    long double deltaT = 1.0E-12;
    /* Скорость выгорания анода */
    long double burningSpeedAnode = 0.44E-5;

    int temperature = 4200;
    int amperage = 150;
    /* Что это */
    double Ti = 4200;
    double Th = 293;
    double Tc = Ti;
    double kb = 1.38E-23;
    double Vemax = sqrt(2 * kb * temperature / massElectron);
    double Vhmax = sqrt(2 * kb * Th / numbersHelium);
    double Vcmax = sqrt(2 * kb * Tc / numbersCarbon);
    double Vhl = Vhmax - 0.1 * Vhmax;
    double Vhpr = Vhmax + 0.1 * Vhmax;
    double Vcl = Vcmax - 0.1 * Vcmax;
    double Vcpr = Vcmax + 0.1 * Vcmax;
    double Vel = Vemax - 0.1 * Vemax;
    double Vepr = Vemax + 0.1 * Vemax;

    /*  Размеры сетки  */
    double xDimensionGrid = 1.0E-3;
    double yDimensionGrid = 1E-2;
    double zDimensionGrid = 1E-2;
    /* число шагов по сетке для крупных частиц */
    int xNumberStepGrid = 5; //im
    int yNumberStepGrid = 5; //km
    int zNumberStepGrid = 5; //lm
    int xInitNumberStepGrid = 10; //imm
    int yInitNumberStepGrid = 10; //kmm
    int zInitNumberStepGrid = 10; //lmm
    /* шаг по сетке */
    double xStepGrid = xDimensionGrid / xNumberStepGrid; //hx
    double yStepGrid = yDimensionGrid / yNumberStepGrid; //hy
    double zStepGrid = zDimensionGrid / zNumberStepGrid; //hz
    /* шаг по сетке(для начального положения частиц) */
    double xInitStepGrid = xDimensionGrid / xInitNumberStepGrid; //hxx
    double yInitStepGrid = yDimensionGrid / yInitNumberStepGrid; //hyy
    double zInitStepGrid = zDimensionGrid / zInitNumberStepGrid; //hzz


    // long double kb = 1.38E-23;

    const int m = xInitNumberStepGrid * yInitNumberStepGrid * zInitNumberStepGrid, n = 3;

    struct Point positionElectron[m][n];
    struct Point positionHelium[m][n];
    struct Point positionCarbon[m][n];
    double speedElectron[m][n];
    double speedHelium[m][n];
    double speedCarbon[m][n];
    
    
    int tt = 0;
    int je = 0;
    for (int i = 1; i <= xInitNumberStepGrid; i++)
    {
        for (int k = 1; k <= yInitNumberStepGrid; k++)
        {
            for (int l = 1; l <= zInitNumberStepGrid; l++)
            {
                // Electron
                // printf("----sdf- %25.20f   ", xInitStepGrid * ((i - 1) + randFrom0To1() / 1E+12)); //xInitNumberStepGrid);
                positionElectron[je][tt].x = xInitStepGrid * ((i - 1) + randFrom0To1() / 1E+12);
                positionElectron[je][tt].y = yInitStepGrid * ((k - 1) + randFrom0To1() / 1E+12);
                positionElectron[je][tt].z = zInitStepGrid * ((l - 1) + randFrom0To1() / 1E+12);
                printf("%35.34f %35.34f %35.34f\n", positionElectron[je][tt].x,  positionElectron[je][tt].y,  positionElectron[je][tt].z );
                // printf("2 %e %e %e\n", xInitStepGrid * ((i - 1) + randFrom0To1() / 1E+12),yInitStepGrid * ((k - 1) + randFrom0To1() / 1E+12),zInitStepGrid * ((l - 1) + randFrom0To1() / 1E+12) );
                speedElectron[je][tt] = sqrt((double)(3  * kb * temperature / massElectron));
                // Helium
                positionHelium[je][tt].x = xInitStepGrid * ((i - 1) + randFrom0To1() / 1E+12);
                positionHelium[je][tt].y = yInitStepGrid * ((k - 1) + randFrom0To1() / 1E+12);
                positionHelium[je][tt].z = zInitStepGrid * ((l - 1) + randFrom0To1() / 1E+12);
                speedHelium[je][tt] = Vhl + (Vhpr - Vhl) * randFrom0To1();
                //Carbon
                positionCarbon[je][tt].x = 0;
                positionCarbon[je][tt].y = yInitStepGrid * ((k - 1) + randFrom0To1() / 1E+12);
                positionCarbon[je][tt].z = zInitStepGrid * ((l - 1) + randFrom0To1() / 1E+12);
                speedCarbon[je][tt] = Vcl + (Vcpr - Vcl) * randFrom0To1();
                // printf("%d\n", je);
                je++;
            }
        }
    }
    // printf("%d %d\n", m, je);
    // struct Point** pe = (struct Point[][]) positionElectron;
    
    // printf("%x %x\n", positionElectron, pe);
    // printf("%x %x\n", positionElectron[0], pe[0]);

    // exit(0);
    for (int i = 1; i <= xNumberStepGrid; i++)
    {
        for (int k = 1; k <= yNumberStepGrid; k++)
        {
            for (int l = 1; l <= zNumberStepGrid; l++)
            {
                // struct Point **pe = (struct Point **)positionElectron;
                for (int j = 0; j < m - 1; j++)
                {
                    printf("%e\n", positionElectron[j][tt].x);
                    if ((pe[j][tt].x >= i * xStepGrid) &&
                            (pe[j][tt].x < (i + 1)*xStepGrid) &&
                            (pe[j][tt].y >= k * yStepGrid) &&
                            (pe[j][tt].y < (k + 1)*yStepGrid) &&
                            (pe[j][tt].z >= l * zStepGrid) &&
                            (pe[j][tt].z < (l + 1)*zStepGrid))
                    {

                    }
                }

            }
        }
    }
    return 0;
}
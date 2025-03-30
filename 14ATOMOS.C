#include <stdio.h>
#include <math.h>

void SuttonChen(float a, int n, int m, float eps, float c, float atm[14][3]);

int main(void)
{
    int n, m;
    float eps, c, a;

    // Parámetros para el níquel
    a = 3.52;
    eps = 0.015707;
    n = 9;
    m = 6;
    c = 39.432;

    float atomos[14][3] = 
    {
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0},
        {1.0, 0.0, 1.0},
        {0.5, 0.0, 0.5},
        {0.0, 1.0, 0.0},
        {1.0, 1.0, 0.0},
        {0.0, 1.0, 1.0},
        {1.0, 1.0, 1.0},
        {0.5, 1.0, 0.5},
        {0.0, 0.5, 0.5},
        {1.0, 0.5, 0.5},
        {0.5, 0.5, 0.0},
        {0.5, 0.5, 1.0}
    };

    SuttonChen(a, n, m, eps, c, atomos);
    return 0;
}

void SuttonChen(float a, int n, int m, float eps, float c, float atm[14][3])
{
    float r[14][14], V[14][14], rho[14], U[14], Utot, x;
    float atm_scaled[14][3];

    for (x = 0.1; x < 10; x += 0.05)
    {
        Utot = 0.0; // Reiniciar energía total en cada iteración

        // Escalar las posiciones de los átomos
        for (int i = 0; i < 14; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                atm_scaled[i][j] = atm[i][j] * x;
            }
        }

        // Inicializar rho y U para cada átomo
        for (int i = 0; i < 14; i++)
        {
            rho[i] = 0.0;
            U[i] = 0.0;
        }

        // Cálculo de distancias, potenciales y densidad electrónica
        for (int i = 0; i < 14; i++)
        {
            for (int j = 0; j < 14; j++)
            {
                if (i != j)
                {
                    r[i][j] = sqrt(pow(atm_scaled[i][0] - atm_scaled[j][0], 2) +
                                   pow(atm_scaled[i][1] - atm_scaled[j][1], 2) +
                                   pow(atm_scaled[i][2] - atm_scaled[j][2], 2));

                    V[i][j] = pow(a / r[i][j], n);
                    rho[i] += pow(a / r[i][j], m);
                }
            }
        }

        // Cálculo de la energía total
        for (int i = 0; i < 14; i++)
        {
            float sum_V = 0.0;
            for (int j = 0; j < 3; j++)
            {
                if (i != j)
                {
                    sum_V += V[i][j];
                }
            }
            U[i] = eps * (0.5 * sum_V - c * sqrt(rho[i]));
            Utot += U[i];
        }

        printf("%f, %f\n", x, Utot);
    }
}

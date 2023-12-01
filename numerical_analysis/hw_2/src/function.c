#include "main.h"

/*
 * task 1: \int _ 0 ^ 1 sin(x) / x dx
 */
double Task1Value(double x)
{
    if (x == 0)
    {
        return 1;
    }
    else
    {
        return sin(x) / x;
    }
}

/*
 * task 2: \int _ 0 ^ 1 e ^ {x ^ 2} dx
 */
double Task2Value(double x)
{
    return exp(x * x);
}

void TrapezoidIntegral(int n, double (*pFun)(double))
{
    double a = 0., b = 1.;
    double h = (b - a) / n;
    double *x = NULL;

    if ((x = (double *)malloc((n + 1) * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! - \'trapezoid method mesh\'\n");
        exit(EXIT_FAILURE);
    }
    for (int index = 0; index <= n; ++index)
    {
        x[index] = a + h * index;
    }

    double value = 0.;
    for (int index = 0; index < n - 1; ++index)
    {
        value += h * (pFun(x[index]) + pFun(x[index + 1])) / 2;
    }
    printf("Trapezoid method: Numerical integration is %021.16le\n", value);

    // free memory
    free(x);
}

void SimpsonIntegral(int n, double (*pFun)(double))
{
    double a = 0., b = 1.;
    double h = (b - a) / n;
    double *x = NULL;

    if ((x = (double *)malloc((n + 1) * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! - \'simpson method mesh\'\n");
        exit(EXIT_FAILURE);
    }
    for (int index = 0; index < n + 1; ++index)
    {
        x[index] = a + h * index;
    }

    double value = 0.;
    for (int index = 0; index < n - 1; ++index)
    {
        value += h * (pFun(x[index]) / 6 +
                      4 * pFun((x[index] + x[index + 1]) / 2) / 6 +
                      pFun(x[index + 1]) / 6);
    }
    printf("Simpson method: Numerical integration is %021.16le\n", value);

    // free memory
    free(x);
}

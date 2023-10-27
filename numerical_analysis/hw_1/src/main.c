#include "main.h"

int main(int argc, char **argv)
{
    puts("==== test ====");
    int n_dim = 0, max_it = 0;
    double r_tol = 0., sor_omega = 0.;
    char *solver_type = NULL;

    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-n_dim", argv[index]))
        {
            n_dim = atoi(argv[index + 1]);
        }
        if (strstr("-max_it", argv[index]))
        {
            max_it = atoi(argv[index + 1]);
        }
        if (strstr("-rtol", argv[index]))
        {
            r_tol = atof(argv[index + 1]);
        }
        if (strstr("-omega", argv[index]))
        {
            sor_omega = atof(argv[index + 1]);
        }
        if (strstr("-solver", argv[index]))
        {
            solver_type = argv[index + 1];
        }
    }
#if 1
    printf("n_dim = %d\n", n_dim);
    printf("max_it = %d, rtol = %.2le\n", max_it, r_tol);
    printf("sor_omega = %.2lf\n", sor_omega);
    printf("solver type = \"%s\"\n", solver_type);
#endif

    double *mat = NULL, *rhs = NULL, *sol = NULL;
    double *ana_sol = NULL; // analytical solution

    LinsysGeneration(n_dim, &mat, &rhs, &sol, &ana_sol);
#if 0
    puts("==== matrix");
    for (int index_i = 0; index_i < n_dim; ++index_i)
    {
        for (int index_j = 0; index_j < n_dim; ++index_j)
        {
            printf("%.4lf\t", mat[index_i * n_dim + index_j]);
        }
        putchar('\n');
    }

    puts("\n==== rhs");
    for (int index = 0; index < n_dim; ++index)
    {
        printf("%.4lf\n", rhs[index]);
    }

    puts("\n==== ana_sol");
    for (int index = 0; index < n_dim; ++index)
    {
        printf("%.4lf\n", ana_sol[index]);
    }

    puts("\n==== sol");
    for (int index = 0; index < n_dim; ++index)
    {
        printf("%.4lf\n", sol[index]);
    }
#endif

    puts("\n==== solver ====");
    if (strstr("lu", solver_type))
    {
        printf("solve with %s\n", solver_type);
        SolverLU(n_dim, mat, rhs, sol);
        printf("|| sol - ana_sol || _ 2 = %021.14le\n", VecNorm(n_dim, sol, ana_sol));
    }
    else if (strstr("jacobi", solver_type))
    {
        printf("solve with %s\n", solver_type);
        SolverJacobi(n_dim, mat, rhs, max_it, r_tol, sol);
        printf("|| sol - ana_sol || _ 2 = %021.14le\n", VecNorm(n_dim, sol, ana_sol));
    }
    else if (strstr("gs", solver_type))
    {
        printf("solve with %s\n", solver_type);
        SolverGS(n_dim, mat, rhs, max_it, r_tol, sol);
        printf("|| sol - ana_sol || _ 2 = %021.14le\n", VecNorm(n_dim, sol, ana_sol));
    }
    else if (strstr("sor", solver_type))
    {
        printf("solve with %s\n", solver_type);
        SolverSOR(n_dim, mat, rhs, max_it, r_tol, sor_omega, sol);
        printf("|| sol - ana_sol || _ 2 = %021.14le\n", VecNorm(n_dim, sol, ana_sol));
    }

    // free memory
    free(mat);
    free(rhs);
    free(sol);
    free(ana_sol);

    return 0;
}

/*
 * command line:
 *     ./app_exe -n <num_1> -max_it <num_2> -rtol <num_3>
 *     -solver <str> -omega <num_4>
 *     <num_1> dimension
 *     <num_2> maximum iteration
 *     <num_3> relative tolerance
 *     <num_4> sor omega
 *     <str>   solver type
 */

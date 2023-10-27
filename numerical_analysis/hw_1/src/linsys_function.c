#include "main.h"

void LinsysGeneration(int dim, double **mat, double **rhs,
                      double **sol, double **ana_sol)
{
    // memory allocation
    /*
     * size mat = (dim x dim, 1)
     * size rhs = (dim, 1)
     * size sol = (dim, 1)
     * size ana_sol = (dim, 1)
     */
    if ((*mat = (double *)malloc(dim * dim * sizeof(double))) == NULL ||
        (*rhs = (double *)malloc(dim * sizeof(double))) == NULL ||
        (*sol = (double *)malloc(dim * sizeof(double))) == NULL ||
        (*ana_sol = (double *)malloc(dim * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed! \"linear system\"\n");
        exit(EXIT_FAILURE);
    }

    // assigning values
    /*
     * ana_sol[i] = 1., i = 1, 2, ..., dim
     * mat[i][j] = 1. / (i + j - 1), i, j = 1, 2, ..., dim
     * rhs = mat * ana_sol
     * sol[i] = 0., i = 1, 2, ..., dim
     */
    for (int index = 0; index < dim; ++index)
    {
        (*ana_sol)[index] = 1.;
    }
    for (int index_i = 0; index_i < dim; ++index_i)
    {
        for (int index_j = 0; index_j < dim; ++index_j)
        {
            (*mat)[index_i * dim + index_j] = 1. / (index_i + index_j + 1);
        }
    }
    for (int index = 0; index < dim; ++index)
    {
        double tmp = 0.;
        for (int index_k = 0; index_k < dim; ++index_k)
        {
            tmp += (*mat)[index * dim + index_k] * (*ana_sol)[index_k];
        }
        (*rhs)[index] = tmp;
    }
    for (int index = 0; index < dim; ++index)
    {
        (*sol)[index] = 0.;
    }
}

void SolverLU(int n, double *mat, double *rhs, double *sol)
{
    // gaussian elimination
    for (int index_i = 0; index_i < n; ++index_i)
    {
        // dominant element in index_i-th column
        for (int index_j = index_i + 1; index_j < n; ++index_j)
        {
            if (mat[index_j * n + index_i] > mat[index_i * n + index_i])
            {
                double tmp = 0.;

                // swap rhs[index_i] rhs[index_j]
                tmp = rhs[index_i];
                rhs[index_i] = rhs[index_j];
                rhs[index_j] = tmp;

                // swap mat[index_i][index_i : end] mat[index_j][index_i : end]
                for (int index_k = index_i; index_k < n; ++index_k)
                {
                    tmp = mat[index_i * n + index_k];
                    mat[index_i * n + index_k] = mat[index_j * n + index_k];
                    mat[index_j * n + index_k] = tmp;
                }
            }
        }

        // elimination
        /*
         * mat[index_j][index_i : end] -= mat[index_j][index_i] / mat[index_i][index_i]
         *                                * mat[index_i][index_i : end]
         * b[index_j] -= mat[index_j][index_i] / mat[index_i][index_i]
         *               * b[index_i]
         */
        for (int index_j = index_i + 1; index_j < n; ++index_j)
        {
            double tmp = mat[index_j * n + index_i] / mat[index_i * n + index_i];
            rhs[index_j] -= tmp * rhs[index_i];
            for (int index_k = index_i; index_k < n; ++index_k)
            {
                mat[index_j * n + index_k] -= tmp * mat[index_i * n + index_k];
            }
        }
    }

    // solving
    /*
     * sol[k] = (rhs[k] - (mat[k][n]*sol[n] + ... mat[k][k + 1] * sol[k + 1])) / mat[k][k],
     *      k = n, ..., 2, 1
     */
    for (int index = n - 1; index >= 0; --index)
    {
        double tmp = 0;
        for (int index_k = n - 1; index_k > index; --index_k)
        {
            tmp += mat[index * n + index_k] * sol[index_k];
        }
        sol[index] = (rhs[index] - tmp) / mat[index * n + index];
    }
}

void SolverJacobi(int n, const double *mat, const double *rhs,
                  int max_it, double r_tol, double *sol)
{
    double *sol_tmp = NULL, *residual = NULL, *vec_zero = NULL;

    // memory allocation
    /*
     * size sol_tmp = n x 1
     * size residual = n x 1
     * size vec_zero = n x 1
     */
    if ((sol_tmp = (double *)malloc(n * sizeof(double))) == NULL ||
        (residual = (double *)malloc(n * sizeof(double))) == NULL ||
        (vec_zero = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"jacobian iteration\"\n");
        exit(EXIT_FAILURE);
    }
    for (int index = 0; index < n; ++index)
    {
        sol_tmp[index] = 0.;  // initialize sol_tmp = 0
        vec_zero[index] = 0.; // all zero vector
    }

    int iter = 0;
    double iter_res = 1.;
    double norm_b = VecNorm(n, rhs, vec_zero);
    printf("\niter\t||true residual||\t||relative residual||\n");
    while (iter <= max_it && iter_res > r_tol)
    {
        // computing residual
        /*
         * r = b - A x
         *     1. computing r = A x
         *     2. computing r = b - r
         */
        MatByVecProduct(n, mat, sol, residual);
        for (int index = 0; index < n; ++index)
        {
            residual[index] = rhs[index] - residual[index];
        }
        double norm_res = VecNorm(n, residual, vec_zero);
        iter_res = norm_res / norm_b;
        printf("%d\t%021.14le\t%021.14le\n", iter, norm_res, iter_res);
        ++iter;

        // jacobian iteration
        /*
         * matrix splitting: A = D - N
         * D * sol = N * sol_tmp + b
         * A[i][i] * sol[i] = -\sum (A[i][j] * sol_tmp[j], j = 1, ..., i - 1)
         *                    -\sum (A[i][j] * sol_tmp[j], i = i + 1, ..., n)
         *                    + b[i]
         */
        JacobiIteration(n, mat, sol_tmp, rhs, sol);

        // updating sol_tmp
        /*
         * sol_tmp[i] = sol[i], i = 1, 2, ..., n
         */
        for (int index = 0; index < n; ++index)
        {
            sol_tmp[index] = sol[index];
        }
    }

    // free memory
    free(sol_tmp);
    free(residual);
    free(vec_zero);
}

void JacobiIteration(int n, const double *mat, const double *x, const double *b, double *y)
{
    for (int index_i = 0; index_i < n; ++index_i)
    {
        double tmp = 0.;
        for (int index_j = 0; index_j < index_i; ++index_j)
        {
            tmp += mat[index_i * n + index_j] * x[index_j];
        }
        for (int index_j = index_i + 1; index_j < n; ++index_j)
        {
            tmp += mat[index_j * n + index_j] * x[index_j];
        }
        tmp = b[index_i] - tmp;
        y[index_i] = tmp / mat[index_i * n + index_i];
    }
}

void MatByVecProduct(int n, const double *mat, const double *x, double *y)
{
    // computing y = mat x
    /*
     * y[index_i] = \sum mat[index_i][index_k] * x[index_k], index_k = 1, ..., n
     */
    for (int index_i = 0; index_i < n; ++index_i)
    {
        double tmp = 0.;
        for (int index_j = 0; index_j < n; ++index_j)
        {
            tmp += mat[index_i * n + index_j] * x[index_j];
        }
        y[index_i] = tmp;
    }
}

void SolverGS(int n, const double *mat, const double *rhs,
              int max_it, double r_tol, double *sol)
{
    double *sol_tmp = NULL, *residual = NULL, *vec_zero = NULL;

    // memory allocation
    /*
     * size sol_tmp = n x 1
     * size residual = n x 1
     * size vec_zero = n x 1
     */
    if ((sol_tmp = (double *)malloc(n * sizeof(double))) == NULL ||
        (residual = (double *)malloc(n * sizeof(double))) == NULL ||
        (vec_zero = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"gauss-seidel iteration\"\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n; ++index)
    {
        sol_tmp[index] = 0.;  // initialize sol_tmp
        vec_zero[index] = 0.; // all zero vector
    }

    int iter = 0;
    double norm_b = VecNorm(n, rhs, vec_zero);
    double iter_res = 1.;
    printf("\niter\t||true residual||\t||relative residual||\n");
    while (iter <= max_it && iter_res > r_tol)
    {
        // computing residual
        /*
         * r = b - A x
         *     1. computing r = A x
         *     2. computing r = b - r
         */
        MatByVecProduct(n, mat, sol, residual);
        for (int index = 0; index < n; ++index)
        {
            residual[index] = rhs[index] - residual[index];
        }
        double norm_res = VecNorm(n, residual, vec_zero);
        iter_res = norm_res / norm_b;
        printf("%d\t%021.14le\t%021.14le\n", iter, norm_res, iter_res);
        ++iter;

        // gs iteration
        /*
         * matrix splitting A = M - N
         *     1. M = D - L, L is lower triangular part of A
         *     2. M * sol = N * sol_tmp + b
         *     3. A[i][i] * sol[i] = -(A[i][j] * sol[j], j = 1, ..., i - 1)
         *                           -(A[i][j] * sol_tmp[j], j = i + 1, ..., n)
         *                           + b[i]
         */
        GSIteration(n, mat, sol_tmp, rhs, sol);

        // updating sol_tmp
        /*
         * sol_tmp[i] = sol[i], i = 1, 2, ..., n
         */
        for (int index = 0; index < n; ++index)
        {
            sol_tmp[index] = sol[index];
        }
    }

    // free memory
    free(sol_tmp);
    free(residual);
    free(vec_zero);
}

void GSIteration(int n, const double *mat, const double *x, const double *b, double *y)
{
    for (int index_i = 0; index_i < n; ++index_i)
    {
        double tmp = 0.;
        for (int index_j = 0; index_j < index_i; ++index_j)
        {
            tmp += mat[index_i * n + index_j] * y[index_j];
        }
        for (int index_j = index_i + 1; index_j < n; ++index_j)
        {
            tmp += mat[index_i * n + index_j] * x[index_j];
        }
        tmp = b[index_i] - tmp;
        y[index_i] = tmp / mat[index_i * n + index_i];
    }
}

void SolverSOR(int n, const double *mat, const double *rhs,
               int max_it, double r_tol, double omega, double *sol)
{
    double *sol_tmp = NULL, *residual = NULL, *vec_zero = NULL;

    // memory allocation
    /*
     * size sol_tmp = n x 1
     * size residual = n x 1
     * size vec_zero = n x 1
     */
    if ((sol_tmp = (double *)malloc(n * sizeof(double))) == NULL ||
        (residual = (double *)malloc(n * sizeof(double))) == NULL ||
        (vec_zero = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"sor iteration\"\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n; ++index)
    {
        sol_tmp[index] = 0.;  // initialize sol_tmp
        vec_zero[index] = 0.; // all zero vector
    }

    int iter = 0;
    double norm_b = VecNorm(n, rhs, vec_zero);
    double iter_res = 1.;
    while (iter <= max_it && iter_res > r_tol)
    {
        // computing residual r = b - A x
        /*
         * 1. r = A x
         * 2. r = b - r
         */
        MatByVecProduct(n, mat, sol, residual);
        for (int index = 0; index < n; ++index)
        {
            residual[index] = rhs[index] - residual[index];
        }
        double norm_res = VecNorm(n, residual, vec_zero);
        iter_res = norm_res / norm_b;
        printf("%d\t%021.14le\t%021.14le\n", iter, norm_res, iter_res);
        ++iter;

        // sor iteration
        /*
         * matrix splitting A = M - N
         *     1. M = 1 / omega D - L
         *     2. (D - omega L) sol = (omega U + (1 - omega) D) sol_tmp + omega b
         *     3. sol[i] = (1 - omega)sol_tmp[i]
         *                 + omega (-(A[i][j] sol[j], j = 1, ..., i - 1)
         *                          -(A[i][j] sol[j], j = i + 1, ..., n)
         *                          + b[i]) / A[i][i]
         */
        SORIteration(n, mat, sol_tmp, rhs, omega, sol);

        // updating sol_tmp
        /*
         * sol_tmp[i] = sol[i], i = 1, 2, ..., n
         */
        for (int index = 0; index < n; ++index)
        {
            sol_tmp[index] = sol[index];
        }
    }

    // free memory
    free(sol_tmp);
    free(residual);
    free(vec_zero);
}

void SORIteration(int n, const double *mat, const double *x, const double *b,
                  double omega, double *y)
{
    for (int index_i = 0; index_i < n; ++index_i)
    {
        double tmp = 0.;
        for (int index_j = 0; index_j < index_i; ++index_j)
        {
            tmp += mat[index_i * n + index_j] * y[index_j];
        }
        for (int index_j = index_i + 1; index_j < n; ++index_j)
        {
            tmp += mat[index_i * n + index_j] * x[index_j];
        }
        tmp = b[index_i] - tmp;
        tmp /= mat[index_i * n + index_i];
        y[index_i] = (1 - omega) * x[index_i] + omega * tmp;
    }
}

double VecNorm(int n, const double *a, const double *b)
{
    double tmp = 0;
    for (int index = 0; index < n; ++index)
    {
        tmp += (a[index] - b[index]) * (a[index] - b[index]);
    }

    return sqrt(tmp);
}

#if 0
void LUDecomposition(int n, const double *mat, double **mat_l, double **mat_u)
{
    // memory allocation
    /*
     * lower triangular matrix mat_l: n x (n + 1) / 2
     * upper triangular matrix mat_u: n x (n + 1) / 2
     */
    if ((*mat_l = (double *)malloc(n * (n + 1) / 2 * sizeof(double))) == NULL ||
        (*mat_u = (double *)malloc(n * (n + 1) / 2 * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"lu decomposition\"\n");
        exit(EXIT_FAILURE);
    }

    // row/col_ptr for matrix mat_l, mat_u
    /*
     * size row_ptr_l = (n + 1) x 1
     * size col_ptr_u = (n + 1) x 1
     */
    int *row_ptr_l = NULL, *col_ptr_u = NULL;
    if ((row_ptr_l = (int *)malloc((n + 1) * sizeof(int))) == NULL ||
        (col_ptr_u = (int *)malloc((n + 1) * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\"row_ptr lu decomposition\"\n");
        exit(EXIT_FAILURE);
    }

    row_ptr_l[0] = 0;
    col_ptr_u[0] = 0;
    for (int index = 1; index <= n; ++index)
    {
        int cnt_l = index; // count of nnz in index-th row
        int cnt_u = index; // count of nnz in index-th row
        row_ptr_l[index] = row_ptr_l[index - 1] + cnt_l;
        col_ptr_u[index] = col_ptr_u[index - 1] + cnt_u;
    }

    // assigning values for mat_l
    /*
     * mat_l[i][i] = 1., i = 1, 2, ..., n
     */
    for (int index = 0; index < n; ++index)
    {
        (*mat_l)[row_ptr_l[index + 1] - 1] = 1.;
    }

    // computing mat_l, mat_u
    for( int index_i = 0; index_i < n; ++index_i )
    {
        for( int index_j = 0; index_j < n; ++index_j )
        {
            int index_start_l = row_ptr_l[index_i];
            int index_end_l = row_ptr_l[index_i + 1];
            int index_start_u = col_ptr_u[index_j];
            int index_end_u = col_ptr_u[index_j + 1];

            // if index_i < index_j, computing mat_l
            /*
            * mat_l[index_i][k] = mat[index_i][k]
            */
        }
    }

    // free memory
    free(row_ptr_l);
    free(col_ptr_u);
}
#endif

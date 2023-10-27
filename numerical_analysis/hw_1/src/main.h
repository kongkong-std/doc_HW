#ifndef MAIN_H_
#define MAIN_H_

// header file
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// function prototype
/*
 * generate linear system
 *    dimension
 *    matrix
 *    right-hand side
 *    numerical solution
 *    analytical solution
 */
void LinsysGeneration(int, double **, double **,
                      double **, double **);

/*
 * lu solver
 *     dimension
 *     matrix
 *     right-hand side
 *     numerical solution
 */
void SolverLU(int, double *, double *, double *);

/*
 * jacobian iteartion solver
 *     dimension
 *     matrix
 *     right-hand side
 *     maximum itertion
 *     relative tolerance
 *     numerical solution
 */
void SolverJacobi(int, const double *, const double *, int, double, double *);

/*
 * gauss seidel iteration solver
 *     dimension
 *     matrix
 *     right-hand side
 *     maximum iteration
 *     relative tolerance
 *     numerical solution
 */
void SolverGS(int, const double *, const double *, int, double, double *);

/*
 * sor iteration solver
 *     dimension
 *     matrix
 *     right-hand side
 *     maximum iteration
 *     relative tolerance
 *     sor omega
 *     numeircal solution
 */
void SolverSOR(int, const double *, const double *, int, double, double, double *);

/*
 * computing L2 norm of vec_a - vec_b
 *     dimension
 *     vector a
 *     vector b
 */
double VecNorm(int, const double *, const double *);

/*
 * computing matrix by vector product y = A x
 *     dimension
 *     matrix A
 *     vector x
 *     vector y
 */
void MatByVecProduct(int, const double *, const double *, double *);

/*
 * jacobian iteration scheme D y = N x + b, A = D - N
 *     dimension
 *     matrix A
 *     vector x
 *     vector b
 *     vector y
 */
void JacobiIteration(int, const double *, const double *, const double *, double *);

/*
 * gauss-seidel iteration scheme M y = N x + b, A = M - N, M = D - L
 *     dimension
 *     matrix A
 *     vector x
 *     vector b
 *     vector y
 */
void GSIteration(int, const double *, const double *, const double *, double *);

/*
 * sor iteration scheme M y = N x + b
 *     dimension
 *     matrix A
 *     vector x
 *     vector b
 *     omega
 *     vector y
 */
void SORIteration(int, const double *, const double *, const double *, double, double *);

#if 0
/*
 * lu decomposition
 *     dimension
 *     matrix
 *     lower triangular matrix
 *     upper triangular matrix
 */
void LUDecomposition(int, const double *, double **, double **);
#endif
#endif

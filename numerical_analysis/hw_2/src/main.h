#ifndef MAIN_H_
#define MAIN_H_

// header file
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// function prototype
void TrapezoidIntegral(int, double (*pFun)(double));
void SimpsonIntegral(int, double (*pFun)(double));
double Task1Value(double);
double Task2Value(double);

#endif
#ifndef HEADER_USER_DEFINED_INTERPOLATION
#define HEADER_USER_DEFINED_INTERPOLATION

#include "myInclude.h"
#include "myMatrix.h"

void LagrangeCoef(double x[], double y[], Matrix* coef);
void LagrangeEval(double x[], Matrix* coef, double* in, double* out);

void NewtonCoef(double x[], double y[], Matrix* coef);
void NewtonEval(double x[], double y[], Matrix* coef, double* in, double* out);



#endif
#pragma once

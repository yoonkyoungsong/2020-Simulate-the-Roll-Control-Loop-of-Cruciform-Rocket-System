#ifndef HEADER_USER_DEFINED_COMPLEX
#define HEADER_USER_DEFINED_COMPLEX

#include "myInclude.h"

typedef struct {
	double real;
	double imge;
}Complex;


Complex adderComplex(Complex A, Complex B);
Complex multiComplex(Complex A, Complex B);
Complex powComplex(Complex A, int n);
double magnitudeComplex(Complex A);

#endif

#pragma once

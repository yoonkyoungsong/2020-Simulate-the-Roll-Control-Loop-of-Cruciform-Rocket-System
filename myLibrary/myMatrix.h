#ifndef HEADER_USER_DEFINED_MATRIX
#define HEADER_USER_DEFINED_MATRIX

#include "myInclude.h"


typedef struct {
	double** at;
	int rows, cols;
}Matrix;


Matrix createMathrix(int rows, int cols);

Matrix createZerosMathrix(int rows, int cols);
Matrix createOnesMathrix(int rows, int cols);
Matrix createIMathrix(int n);

void freeMatrix(Matrix matrix);

void printMatrix(char name[], Matrix matrix);
void saveMathrix(char Name[], Matrix matrix);


Matrix adderMatrix(Matrix A, Matrix B);
Matrix subtrMatrix(Matrix A, Matrix B);
Matrix multiMatrix(Matrix A, Matrix B);

Matrix multiConstantandMatrix(double k, Matrix A);

double dotProductMatrix(Matrix A, Matrix B);
double magMatrix(Matrix A);
double distanceBiMatrix(Matrix A, Matrix B);

void invers3by3Matrix(Matrix* A);
void det3byMatrix(Matrix* A, double* det);

float determinantMatrix(float a[25][25], float k);
void cofactorMatrix(float num[25][25], float f, float inverse[25][25]);
void transposeMatrix(float num[25][25], float fac[25][25], float r, float inverse[25][25]);

Matrix inverseMatrix(Matrix A);
Matrix transMatrix(Matrix A);

#endif

#pragma once

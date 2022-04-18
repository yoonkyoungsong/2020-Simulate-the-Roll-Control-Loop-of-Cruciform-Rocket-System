#include "myInterpolation.h"


//---------------------Lagrange Interpolation--------------------//
void LagrangeCoef(double x[], double y[], Matrix* coef) {

	double d;

	for (int i = 0; i < (*coef).rows; i++) { (*coef).at[i][0] = 0.0; }

	for (int k = 0; k < (*coef).rows; k++)
	{	
		d = 1;
		
		for (int i = 0; i < (*coef).rows; i++)
		{
			if (i != k) {	d = d * (x[k] - x[i]);	}
			(*coef).at[k][0] = y[k] / d;
		}
	}

}

void LagrangeEval(double x[], Matrix* coef, double* in, double* out) {
	
	double N = 0 ;
	
	for (int k = 0; k < (*coef).rows; k++)
	{
		N = 1;

		for (int i = 0; i < (*coef).rows; i++)
		{
			if (i != k) { N = N * (*in - x[i]); }	
		}

		*out = *out + N * (*coef).at[k][0];
	}
}



//---------------------Newton Interpolation--------------------//

void NewtonCoef(double x[], double y[], Matrix* coef) {

	for (int i = 0; i < (*coef).rows; i++) { (*coef).at[i][0] = y[i]; }

	for (int j = 1; j < (*coef).rows; j++)
	{
		for (int i = 0; i < (*coef).rows - j; i++)
		{
			(*coef).at[i][j] = ( (*coef).at[i+1][j-1] - (*coef).at[i][j-1] ) / (x[i+j] - x[i]);
		}
	}
}

void NewtonEval(double x[], double y[], Matrix* coef, double* in, double* out) {

	double buf = 1;

	NewtonCoef(x, y, coef);

	*out = (*coef).at[0][0];

	for (int i = 1; i < (*coef).rows; i++)
	{
		buf = buf * (*in - x[i - 1]);
		*out = *out + (*coef).at[0][i] * buf;
	}

}


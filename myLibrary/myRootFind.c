#include "myRootFind.h"

unsigned int method;
unsigned int iter = 0;

//----------Bisection variable----------//
double		 a, b, m, bound;

//-----------Newton variable------------//
double		 F[2] = { 0 };
double		 x[MAX_ITER] = {0};

//-----------Secant variable------------//
double		 y[MAX_ITER] = { 0 };
double		 Dx[MAX_ITER] = { 0 };


void EvalFx(double x, double ValFx[]) {

	ValFx[FX] = pow(x, 3.0) - 3.0 * x + 1.0;
	ValFx[DFX] = 3.0 * pow(x, 2.0) - 3.0;

}

int ReadMethod() {

	printf("This is program to find root \n\n");
	printf("1: Bisection   2:Newton   3:Secant \n\n");
	printf("Enter desired root find method number:");
	
	scanf("%d", &method);
	
	printf("\n");

	return method;
}

void RootFinding(int method, int max_iter, double *sol, int *flagError) {

	switch (method) {
	
	case BISECT:
		printf("Bisection Start\n");
		Bisection(max_iter, sol, flagError);
		break;

	case NEWTON:
		printf("Newton Start\n");
		Newton(max_iter, sol, flagError);
		break;

	case SECANT:
		printf("Secant Start\n");
		Secant(max_iter, sol, flagError);
		break;

	default: 
		printf("Another number has been entered\n");
	
	}

}

void Bisection(int max_iter, double* sol, int* flagError) {
	
	a = MIN_ROOT;
	b = MAX_ROOT;

	if(sign(F(a)) == sign(F(b))){ printf("function has same sign at end point\n");}
	else{
		printf("-------------------------------------------------------------\n");
		printf("  step    a         b         m         ym         bound \n");
		
		for (int k = 1; k <= MAX_ITER; k++)
		{	
			m = (a + b) / 2;
			bound = (b - a) / 2;
			iter = k;

			printf("   %d     %.4f    %.4f    %.4f    %.4f    %.4f \n", iter, a, b, m, F(m), bound);

			if (F(m) * F(a) < 0) { b = m; }
			else { a = m; }
			
			if (fabs(F(m)) < TOL_ERROR) 
			{
				*flagError = CONVERGE;
				break;
			}

			if (iter == max_iter)
			{
				*flagError = ERR_ITER;
				break;
			}
		}

		*sol = m;
		printf("-------------------------------------------------------------\n");
	}
}

void Newton(int max_iter, double* sol, int* flagError) {

	iter = 1;
	x[0] = INIT_X;
	EvalFx(x[0], F);

	printf("-------------------------------------------------------------\n");
	printf("  step    x        f(x) \n");
	printf("   %d     %.4f    %.4f \n", iter, x[0], F[FX]);

	if (F[DFX] == 0) { *flagError = ERR_DIVERGE; }
	else
	{

		for (int k = 1; k < MAX_ITER; k++)
		{
			x[k] = x[k - 1] - F[FX] / F[DFX];
			EvalFx(x[k], F);

			printf("   %d     %.4f    %.4f \n", iter + 1, x[k], F[FX]);

			if (fabs(x[k] - x[k - 1]) < TOL_ERROR)
			{
				*flagError = CONVERGE;
				break;
			}
			if (iter == max_iter)
			{
				*flagError = ERR_ITER;
				break;
			}
			if (F[DFX] == 0) {
				*flagError = ERR_DIVERGE;
				break;
			}
			iter++;

		}
		*sol = x[iter - 1];
	}
	printf("-------------------------------------------------------------\n");
}

void Secant(int max_iter, double* sol, int* flagError) {
	
	iter = 1;

	x[0] = INIT_X1;
	x[1] = INIT_X2;
	y[0] = F(x[0]);
	y[1] = F(x[1]);

	printf("-------------------------------------------------------------\n");
	printf("  step    x(k-1)     x(k)      x(k+1)     y(k+1)      Dx(k+1) \n");

	for (int k = 1; k < MAX_ITER; k++)
	{
		x[k + 1] = x[k] - y[k] * (x[k] - x[k - 1]) / (y[k] - y[k - 1]);
		y[k + 1] = F(x[k + 1]);
		Dx[k + 1] = x[k + 1] - x[k];
		iter = k;
		
		printf("   %d     %.4f     %.4f     %.4f     %.4f    %.4f \n", iter, x[k-1], x[k], x[k+1], y[k+1], Dx[k+1]);

		if (fabs(y[k+1]) < TOL_ERROR)
		{
			*flagError = CONVERGE;
			break;
		}
		if (iter == max_iter)
		{
			*flagError = ERR_ITER;
			break;
		}

	}
	*sol = x[iter+1];
	printf("-------------------------------------------------------------\n");

}

void ViewResultRootFinding(double sol) {

	printf("\nRoot finding solution     x = %.4f\n\n", sol);

}

void ErrorMessage(int flagError){

	switch (flagError) {
	
	case CONVERGE:
		printf("Complete: root finding has converged \n");
		break;
		
	case ERR_ITER:
		printf("Error message: zero not found to desired toleranve\n");
		break;

	case ERR_DIVERGE:
		printf("Error message: Zero not found since solution being divergence\n");
		printf("               Reset inital x value\n");
		break;
	
	}


}



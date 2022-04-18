#ifndef HEADER_USER_DEFINED_ROOTFIND
#define HEADER_USER_DEFINED_ROOTFIND

#include "myInclude.h"
#include "myComplex.h"
#include "myMatrix.h"

#define F(x) (pow(x, 3.0) - 3.0 * x + 1.0)
#define sign(x) ((x>0)? 1:-1)

//---------Root finding defination--------//
#define MAX_ITER		(unsigned int)	( 100 )
#define TOL_RESIDUAL	(double)		( 1e-5 )
#define TOL_ERROR		(double)		( 1e-5 )

//----------Bisection defination----------//
#define MIN_ROOT		(double)		( 0.0 )
#define MAX_ROOT		(double)		( 1.0 ) 

//------------Newton defination-----------//
#define INIT_X			(double)		( 0.0 ) 

//------------Secant defination-----------//
#define INIT_X1			(double)		( 0.0 ) 
#define INIT_X2			(double)		( 1.0 ) 


enum RootFindingMethod	{ BISECT=1,		 NEWTON,	SECANT		};
enum FunctionEval		{ FX = 0,		 DFX					};
enum ErrorCode			{ CONVERGE = 0,	 ERR_ITER,  ERR_DIVERGE };
enum ViewRootFinding	{ ONLYSOL = 0,   ALLVIEW				};

void EvalFx(double x, double ValFx[]);

int ReadMethod();

void RootFinding(int method, int max_iter, double* sol, int* flagError);

void Bisection	(int max_iter, double* sol, int* flagError);
void Newton		(int max_iter, double* sol, int* flagError);
void Secant		(int max_iter, double* sol, int* flagError);

void ViewResultRootFinding	(double sol);
void ErrorMessage			(int flagError);

void calculateFunction(Matrix* num, double wgc, double ValFx[]);


#endif

#pragma once

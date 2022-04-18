/*-----------------------------------------------
		   2020-1  Numerucal Analysis
			   Final: ODE SOLVER
		   21600372  Yoonkyoung Song
-----------------------------------------------*/

#include "myInclude.h"
#include "myUserDefInclude.h"
#include "NAFinal.h"


int sysnum = 8;
int iter = 0;

Matrix sys1A, sys1B;
Matrix sys2A, sys2B;
Matrix dstate, state;
Matrix beforeState, beforedstate;
Matrix k;
Matrix delC;

double t[DATANUM] = { 0.0 };
double phic[DATANUM] = {  INITpic };
double delc[DATANUM] = { 0.0 };
double del[DATANUM] = { 0.0 };
double p[DATANUM] = { 0.0 };
double phi[DATANUM] = { 0.0 };


void main()
{
	initalizeSystme(&sys1A, &sys1B, &sys2A, &sys2B);
	initalizeState(&dstate, &state, &beforeState, &beforedstate, &delC, &k, sysnum);
	
	do {
		
		input(phic,iter);

		
		system4(&sys1A, &sys1B, &dstate, &state, &beforeState, &beforedstate);
		system1(&dstate, &state, &beforeState, phic, &delC, iter);
		system2(&sys2A, &sys2B, &dstate, &state, &delC);
		systme3(&dstate, &state, &beforeState);
	

		storageData(&state, &beforeState, &dstate, &beforedstate, &delC, delc, del, p, phi, iter);
		
		integration(&dstate, &state, &k, RK4, &iter);
		checkStop(t, iter);

	} while (t[iter]<tf);

	WriteTxtFile("in0_wa24.txt", DATANUM, t, phic, delc, del, p, phi);
	
}
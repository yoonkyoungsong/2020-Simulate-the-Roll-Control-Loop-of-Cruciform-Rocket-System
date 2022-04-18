#ifndef HEADER_USER_DEFINED_NAFINAL
#define HEADER_USER_DEFINED_NAFINAL

#include "myInclude.h"
#include "myMatrix.h"

#define  PI (double) 3.14159265358979323846

#define INITpic (double) 0.0  // PI/8 PI/4 0

#define Za (double) 0.7
#define Wa (double) 2*PI*24.6 //16.4 24.6 20.5

#define Zs (double) 0.7
#define Ws (double) 2*PI*100

#define Kp (double) 27.0
#define Kd (double) 0.003

#define Lp (double) 3.2
#define Lpi (double) 1200
#define Ldel (double) 16000

#define ts (double) 0.0005
#define tf (double) 2.0

#define DATANUM (int) 4000

enum methodODE { EULER = 1, RK4};
//enum systemVari { DEL = 0, dDEL, PHI, dPHI, hPHI, dhPHI, hP, dhP};

void initalizeState(Matrix* dstate, Matrix* state, Matrix* beforeState, Matrix* beforedstate, Matrix* dleC, Matrix* k, int sysnum);
void initalizeSystme(Matrix* sys1A, Matrix* sys1B, Matrix* sys2A, Matrix* sys2B);

void input(double pic[], int iter);

void system1(Matrix* dstate, Matrix* state, Matrix* beforeState, double input[], Matrix* delC, int iter);
void system2(Matrix* A, Matrix* B, Matrix* dstate, Matrix* state, Matrix* delC);
void systme3(Matrix* dstate, Matrix* state, Matrix* beforeState);
void system4(Matrix* A, Matrix* B, Matrix* dstate, Matrix* state, Matrix* beforeState, Matrix* beforedstate);

void storageData(Matrix* state, Matrix* beforeState, Matrix* dstate, Matrix* beforedstate, Matrix* delC, double delc[], double del[], double p[], double phi[], int iter);

void integration(Matrix* dstate, Matrix* state, Matrix* k, int method, int* iter);

void Euler(Matrix* dstate, Matrix* state);
void RungeKutta4(Matrix* dstate, Matrix* state, Matrix* k);


void checkStop(double t[], int iter);
void WriteTxtFile(char Name[], double N_DATA, double data1[], double data2[], double data3[], double data4[], double data5[], double data6[]);


#endif
#pragma once

#include "NAFinal.h"

void initalizeState(Matrix* dstate, Matrix* state, Matrix* beforeState, Matrix* beforedstate, Matrix* dleC, Matrix* k, int sysnum) {

	*dstate = createZerosMathrix(4, sysnum);
	*state = createZerosMathrix(4, sysnum);
	*beforeState = createZerosMathrix(4, sysnum);
	*beforedstate = createZerosMathrix(4, sysnum);
	*k = createZerosMathrix(4, sysnum);

	*dleC = createZerosMathrix(4, 1);

	state->at[0][2] = INITpic + 1.0 * PI / 180;
}

void initalizeSystme(Matrix* sys1A, Matrix* sys1B, Matrix* sys2A, Matrix* sys2B) {

	*sys1A = createZerosMathrix(2, 2);
	*sys1B = createZerosMathrix(2, 1);
	
	*sys2A = createZerosMathrix(2, 2);
	*sys2B = createZerosMathrix(2, 1);

	sys1A->at[0][0] = 0;
	sys1A->at[0][1] = 1;
	sys1A->at[1][0] = pow(Ws, 2) * (-1);
	sys1A->at[1][1] = Zs * Ws * (-2);

	sys1B->at[0][0] = 0;
	sys1B->at[1][0] = pow(Ws, 2);

	sys2A->at[0][0] = 0;
	sys2A->at[0][1] = 1;
	sys2A->at[1][0] = pow(Wa, 2) * (-1);
	sys2A->at[1][1] = Za * Wa * (-2);

	sys2B->at[0][0] = 0;
	sys2B->at[1][0] = pow(Wa, 2);
}

void input(double pic[], int iter) {

	pic[iter] = INITpic;

}




void system1(Matrix* dstate, Matrix* state, Matrix* beforeState, double input[], Matrix* delC, int iter) {

	for (int i = 0; i < 4; i++) {
		delC->at[i][0] = Kd * (Kp * (input[iter] - state->at[i][6]) - state->at[i][4]);
	}

}

void system2(Matrix* A, Matrix* B, Matrix* dstate, Matrix* state, Matrix* delC) {

	for (int i = 0; i < 4; i++) {

		dstate->at[i][0] = state->at[i][1];
		dstate->at[i][1] = A->at[1][0] * state->at[i][0] + A->at[1][1] * state->at[i][1] + B->at[1][0] * delC->at[i][0];
	}

}

void systme3(Matrix* dstate, Matrix* state, Matrix* beforeState) {

	for (int i = 0; i < 4; i++) {

		dstate->at[i][2] = state->at[i][3];
		dstate->at[i][3] = Lpi * sin(4 * beforeState->at[i][2]) - Lp * state->at[i][3] + Ldel * state->at[i][0];


	}

}

void system4(Matrix* A, Matrix* B, Matrix* dstate, Matrix* state, Matrix* beforeState, Matrix* beforedstate) {


	for (int i = 0; i < 4; i++) {

		dstate->at[i][4] = state->at[i][5];
		dstate->at[i][5] = A->at[1][0] * state->at[i][4] + A->at[1][1] * state->at[i][5] + B->at[1][0] * state->at[i][3];

		dstate->at[i][6] = state->at[i][7];
		dstate->at[i][7] = A->at[1][0] * state->at[i][6] + A->at[1][1] * state->at[i][7] + B->at[1][0] * state->at[i][2];

	}

}



void integration(Matrix* dstate, Matrix* state, Matrix* k, int method, int* iter) {

	switch(method){
	
	case 1:
		Euler(dstate, state);
	
	case 2:

		RungeKutta4(dstate, state, k);

	default :
		break;

	}

	*iter += 1;

}

void Euler(Matrix* dstate, Matrix* state) {

	for (int i = 0; i < (*dstate).cols; i++) {

		state->at[0][i] = state->at[0][i] + ts * dstate->at[0][i];
	}

}

void RungeKutta4(Matrix* dstate, Matrix* state, Matrix* k) {

	for (int j = 0; j < (*dstate).cols/2; j++) {
		
		for (int i = 0; i < 4; i++) {
			k->at[i][2 * j] = ts * dstate->at[i][2 * j + 1];
		}

		dstate->at[0][2 * j] = dstate->at[0][2 * j] + (k->at[0][2 * j] + 2 * k->at[1][2 * j] + 2 * k->at[2][2 * j] + k->at[3][2 * j]) / 6;
		dstate->at[1][2 * j] = dstate->at[0][2 * j] + 0.5 * k->at[0][2 * j];
		dstate->at[2][2 * j] = dstate->at[0][2 * j] + 0.5 * k->at[1][2 * j];
		dstate->at[3][2 * j] = dstate->at[0][2 * j] + k->at[2][2 * j];

		state->at[0][2 * j + 1] = dstate->at[0][2 * j];
		state->at[1][2 * j + 1] = dstate->at[1][2 * j];
		state->at[2][2 * j + 1] = dstate->at[2][2 * j];
		state->at[3][2 * j + 1] = dstate->at[3][2 * j];

	}

	for (int j = 0; j < (*dstate).cols / 2; j++) {

		for (int i = 0; i < 4; i++) {
			k->at[i][2 * j] = ts * dstate->at[i][2 * j];
		}

		state->at[0][2 * j] = state->at[0][2 * j] + (k->at[0][2 * j] + 2 * k->at[1][2 * j] + 2 * k->at[2][2 * j] + k->at[3][2 * j]) / 6;
		state->at[1][2 * j] = state->at[0][2 * j] + 0.5 * k->at[0][2 * j];
		state->at[2][2 * j] = state->at[0][2 * j] + 0.5 * k->at[1][2 * j];
		state->at[3][2 * j] = state->at[0][2 * j] + k->at[2][2 * j];

	}




}

void checkStop(double t[], int iter) {

	t[iter] = t[iter-1] + ts;

}

void storageData(Matrix* state, Matrix* beforeState, Matrix* dstate, Matrix* beforedstate, Matrix* delC, double delc[], double del[], double p[], double phi[], int iter) {

	del[iter] = state->at[0][0];
	delc[iter] = delC->at[0][0];
	p[iter] = state -> at[0][3];
	phi[iter] = state->at[0][2];
	
	*beforeState = *state;
	*beforedstate = *dstate;

}



void WriteTxtFile(char Name[], double N_DATA, double data1[], double data2[], double data3[], double data4[], double data5[], double data6[]) {

	FILE* pF1; 

	//char way[1000] = "C:\\Users\\xde12\\OneDrive\\πŸ≈¡ »≠∏È\\";


	pF1 = fopen(Name, "w");

	for (int count = 0; count < N_DATA; count++)
	{
		fprintf(pF1, "%lf %lf %lf %lf %lf %lf \n", data1[count], data2[count], data3[count], data4[count], data5[count], data6[count]);
	}

	fclose(pF1);

	printf("Storage data file");
}
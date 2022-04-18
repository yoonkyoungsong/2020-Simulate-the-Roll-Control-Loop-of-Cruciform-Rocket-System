#include "myMatrix.h"


Matrix createMathrix(int rows, int cols) {

	Matrix matrix;
	matrix.rows = rows;
	matrix.cols = cols;
	
	matrix.at = malloc(sizeof(double*)*rows);

	for (int i = 0; i < rows; i++)
	{
		matrix.at[i] = malloc(sizeof(double) *cols);
	}

	return matrix;

}

Matrix createZerosMathrix(int rows, int cols) {
	
	Matrix matrix;
	
	matrix = createMathrix(rows, cols);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			matrix.at[i][j] = 0;
		}
	}

	return matrix;
}

Matrix createOnesMathrix(int rows, int cols) {

	Matrix matrix;

	matrix = createMathrix(rows, cols);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			matrix.at[i][j] = 1;
		}
	}

	return matrix;

}

Matrix createIMathrix(int n) {

	Matrix I;
	
	I = createZerosMathrix(n,n);

	for (int i = 0; i < n; i++) { I.at[i][i] = 1; }

	return I;

}

void freeMatrix(Matrix matrix)
{
	for (int i = 0; i < matrix.rows; i++)
	{
		free((matrix.at[i]));
	}

	free(matrix.at);

}

void printMatrix(char name[], Matrix matrix) {

	printf("\n\n");
	printf(name);
	printf(" =\n");
	
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{
			printf( "%f", matrix.at[i][j]);
			printf("     ");
		}
		printf("\n");
	}
	printf("\n\n");
}

void saveMathrix(char Name[], Matrix matrix) {

	FILE* pF1;
	char s1[50];
	char s2[10] = ".txt";

	strcpy(s1, Name);
	strcat(s1, s2);

	pF1 = fopen(Name, "w");

	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{
			fprintf(pF1, "%f ", matrix.at[i][j]);
		}
		fprintf(pF1, "\n");
	}

	fclose(pF1);

	printf("Save matrix file\n");

}

Matrix adderMatrix(Matrix A, Matrix B) {

	Matrix add;

	add = createZerosMathrix(A.rows, A.cols);

	if (A.rows == B.rows && A.cols == B.cols) {

		for (int i = 0; i < add.rows; i++) {
			for (int j = 0; j < add.cols; j++) {
				add.at[i][j] = A.at[i][j] + B.at[i][j];
			}
		}

	}
	else {
		printf("Add:The dimensions of the two matrices that you compute do not match\n\n");
	}

	return add;

}

Matrix subtrMatrix(Matrix A, Matrix B) {

	Matrix sub;

	sub = createZerosMathrix(A.rows, A.cols);

	if (A.rows == B.rows && A.cols == B.cols) {

		for (int i = 0; i < sub.rows; i++) {
			for (int j = 0; j < sub.cols; j++) {
				sub.at[i][j] = A.at[i][j] - B.at[i][j];
			}
		}

	}
	else {
		printf("Sub: The dimensions of the two matrices that you compute do not match\n\n");
	}

	return sub;

}

Matrix multiMatrix(Matrix A, Matrix B) {

	Matrix multi;

 	multi = createZerosMathrix(A.rows, B.cols);

	if (A.cols == B.rows) {

		for (int i = 0; i < multi.rows; i++) {
			for (int j = 0; j < multi.cols; j++) {
				for (int k = 0; k < A.cols; k++) {
					multi.at[i][j] = multi.at[i][j] + (A.at[i][k] * B.at[k][j]);
				}
			}
		}

	}
	else {
		printf("Multi: The dimensions of the two matrices that you compute do not match\n\n");
	}

	return multi;

}

Matrix multiConstantandMatrix(double k, Matrix A) {

	Matrix multi;

	multi = createZerosMathrix(A.rows, A.cols);


	for (int i = 0; i < multi.rows; i++) {
		for (int j = 0; j < multi.cols; j++) {
			multi.at[i][j] = k * A.at[i][j];
		}
	}

	return multi;
}

double magMatrix(Matrix A) {

	double mag = 0;

	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			mag = mag + pow(A.at[i][j], 2);
		}
	}
	return mag = sqrt(mag);
}

double dotProductMatrix(Matrix A, Matrix B) {

	double dot = 0;

	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			dot = A.at[i][j] + B.at[i][j];
		}
	}
	return dot;

}

double distanceBiMatrix(Matrix A, Matrix B) {
	
	return magMatrix(subtrMatrix(A,B));

}

void invers3by3Matrix(Matrix* A) {

	double det;

	det = (*A).at[0][0] * (*A).at[1][1] * (*A).at[2][2]
		- (*A).at[0][0] * (*A).at[1][2] * (*A).at[2][1]
		- (*A).at[0][1] * (*A).at[1][0] * (*A).at[2][2]
		+ (*A).at[0][1] * (*A).at[1][2] * (*A).at[2][0]
		+ (*A).at[0][2] * (*A).at[1][0] * (*A).at[2][1]
		- (*A).at[0][2] * (*A).at[1][1] * (*A).at[2][0];

	(*A).at[0][0] = ((*A).at[1][1] * (*A).at[2][2] - (*A).at[1][2] * (*A).at[2][1]) / det;
	(*A).at[0][1] = ((*A).at[0][2] * (*A).at[2][1] - (*A).at[0][1] * (*A).at[2][2]) / det;
	(*A).at[0][2] = ((*A).at[0][1] * (*A).at[1][2] - (*A).at[0][0] * (*A).at[1][2]) / det;
	
	(*A).at[1][0] = ((*A).at[1][2] * (*A).at[2][0] - (*A).at[1][0] * (*A).at[2][2]) / det;
	(*A).at[1][1] = ((*A).at[0][0] * (*A).at[2][2] - (*A).at[0][2] * (*A).at[2][0]) / det;
	(*A).at[1][2] = ((*A).at[0][2] * (*A).at[1][0] - (*A).at[0][0] * (*A).at[1][2]) / det;
	
	(*A).at[2][0] = ((*A).at[1][0] * (*A).at[2][1] - (*A).at[1][1] * (*A).at[2][0]) / det;
	(*A).at[2][1] = ((*A).at[0][1] * (*A).at[2][0] - (*A).at[0][0] * (*A).at[2][1]) / det;
	(*A).at[2][2] = ((*A).at[0][0] * (*A).at[1][1] - (*A).at[0][1] * (*A).at[1][0]) / det;

}

void det3byMatrix(Matrix* A, double* det) {

	(*det) = (*A).at[0][0] * (*A).at[1][1] * (*A).at[2][2]
		   - (*A).at[0][0] * (*A).at[1][2] * (*A).at[2][1]
		   - (*A).at[0][1] * (*A).at[1][0] * (*A).at[2][2]
		   + (*A).at[0][1] * (*A).at[1][2] * (*A).at[2][0]
		   + (*A).at[0][2] * (*A).at[1][0] * (*A).at[2][1]
		   - (*A).at[0][2] * (*A).at[1][1] * (*A).at[2][0];

}



float determinantMatrix(float a[25][25], float k)
{
	float s = 1, det = 0, b[25][25];
	int m, n;
	if (k == 1)
	{
		return (a[0][0]);
	}
	else
	{
		det = 0;
		for (int c = 0; c < k; c++)
		{
			m = 0;
			n = 0;
			for (int i = 0; i < k; i++)
			{
				for (int j = 0; j < k; j++)
				{
					b[i][j] = 0;
					if (i != 0 && j != c)
					{
						b[m][n] = a[i][j];
						if (n < (k - 2))
							n++;
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (a[0][c] * determinantMatrix(b, k - 1));
			s = -1 * s;
		}
	}

	return (det);
}

void cofactorMatrix(float num[25][25], float f, float inverse[25][25])
{
	float b[25][25], fac[25][25];
	int m, n;
	for (int q = 0; q < f; q++)
	{
		for (int p = 0; p < f; p++)
		{
			m = 0;
			n = 0;
			for (int i = 0; i < f; i++)
			{
				for (int j = 0; j < f; j++)
				{
					if (i != q && j != p)
					{
						b[m][n] = num[i][j];
						if (n < (f - 2))
							n++;
						else
						{
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = pow(-1, q + p) * determinantMatrix(b, f - 1);
		}
	}
	transposeMatrix(num, fac, f, inverse);
}

/*Finding transpose of matrix*/
void transposeMatrix(float num[25][25], float fac[25][25], float r, float inverse[25][25])
{
	float b[25][25], d;

	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < r; j++)
		{
			b[i][j] = fac[j][i];
		}
	}

	d = determinantMatrix(num, r);

	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < r; j++)
		{
			inverse[i][j] = b[i][j] / d;
		}
	}

	/*
	printf("\nThe inverse of matrix is : \n");

	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < r; j++)
		{
			printf("\t%f", inverse[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	*/
}

Matrix inverseMatrix(Matrix A) {
	
	Matrix inv;
	float buff[25][25], inverse[25][25], d;
	
	inv = createZerosMathrix(A.rows, A.cols);

	for (int i = 0; i < A.rows ; i++) 
	{
		for (int j = 0; j < A.cols; j++)
		{
			buff[i][j] = A.at[i][j];
		}
	
	}
	
	d = determinantMatrix(buff, A.rows);
	
	if (d == 0)
		printf("\nInverse of Entered Matrix is not possible\n");
	else
		cofactorMatrix(buff, A.rows,inverse);

	for (int i = 0; i < A.rows; i++)
	{
		for (int j = 0; j < A.cols; j++)
		{
			inv.at[i][j] = inverse[i][j];
		}

	}
	
	return inv;

}

Matrix transMatrix(Matrix A) {

	Matrix trans;

	trans = createZerosMathrix(A.cols, A.rows);

	for (int i = 0; i < trans.rows; i++) {
		for (int j = 0; j < trans.cols; j++) {
			trans.at[i][j] = A.at[j][i];
		}
	}

	return trans;

}

#include "myApproximation.h"

void LSApproximation_Batch(Matrix* data, int n_data, int order, Matrix* coeff) {
	
	Matrix H,y;

	H = createZerosMathrix(n_data, order);
	y = createZerosMathrix(n_data, 1);

	for (int i = 0; i < H.rows; i++) {	
		for (int j = 0; j < H.cols; j++) 
		{	
			if (j == H.cols - 1) {	H.at[i][j] = 1;	 }
			else {
				H.at[i][j] = pow((*data).at[i][0], H.cols - j - 1);
			}
		}
		y.at[i][0] = (*data).at[i][1];
	}

	(*coeff) = multiMatrix(inverseMatrix(multiMatrix(transMatrix(H), H)), multiMatrix(transMatrix(H), y));

	printf("-------------------- LS Solution ------------------------\n");
	for (int i = 0; i < order; i++) {
		printf("coeff a[%d] = %lf \n", order-i-1,(*coeff).at[i][0]);
	}
	printf("---------------------------------------------------------\n");

}



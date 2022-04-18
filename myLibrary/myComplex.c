#include "myComplex.h"

Complex adderComplex(Complex A, Complex B) {

	Complex add;

	add.real = A.real + B.real;
	add.imge = A.imge + B.imge;

	return add;

}

Complex multiComplex(Complex A, Complex B) {

	Complex multi;

	multi.real = A.real * B.real - A.imge * B.imge ;
	multi.imge = A.imge * B.real + A.real * B.imge;

	return multi;
}

Complex powComplex(Complex A, int n) {

	Complex pow;

	pow = A;

	if (n == 0) { pow.real = 1; pow.imge = 0; }
	else if (n == 1) { }
	else if(n > 1) {
		for (int i = 0; i < n - 1; i++) {
			pow = multiComplex(pow, A);
		}
	}

	return pow;
}

double magnitudeComplex(Complex A) {

	double mag;

	mag = sqrt(pow(A.real,2)+pow(A.imge, 2));

	return mag;

}

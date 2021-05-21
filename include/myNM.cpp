/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"
#include "myMatrix.h"
#include "math.h"

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

// Matrix subtraction
Matrix	subMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] - _B.at[i][j];

	return Out;
}

// Apply forward-substitution
void fwdsub(Matrix _L, Matrix _y, Matrix _b)// Lower triangle과  vector b를 입력받으면 y를 계산한다
{

	for (int i = 0; i < _L.rows ; i++) {
		double temp = 0;

		for (int j = 0 ; j <= i - 1 ; j++)
			temp += _L.at[i][j] * _y.at[j][0];

		_y.at[i][0] = (_b.at[i][0] - temp) / _L.at[i][i];
	}
}

// Apply back-substitution
void backSub(Matrix _U, Matrix _y, Matrix _x)// upper triangle과 변경된 vector d를 입력받으면 x를 계산한다
{

	for (int i = _U.rows-1; i >= 0; i--) {
		
		double temp = 0;
		for (int j = i ; j < _U.cols; j++)
			temp += _U.at[i][j] * _x.at[j][0];
		_x.at[i][0] = (_y.at[i][0] - temp) / _U.at[i][i];
	}
}

void gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d)//입력받은 행렬을 upper triangle 행렬과 변경된 vector d로 변환시키는 함수
{
	for (int i = 0; i < _A.cols; i++) {
		_d.at[i][0] = _b.at[i][0];//입력받은 vector b를 vector d 에 복사
		for (int j = 0; j < _A.rows; j++) {
			_U.at[i][j] = _A.at[i][j];//입력받은 Matrix U를 Matrix A 에 복사
		}
	}

	for (int k = 0; k < _A.cols-1; k++) {
		for (int i = k + 1; i < _A.cols; i++) {
				_d.at[i][0] = _d.at[i][0] - (_U.at[i][k] / _U.at[k][k]) * _d.at[k][0];
				double temp = _U.at[i][k];// 아래 반복문에서 계산이 처음 이뤄지는 항이 0이되면 나머지 계산을 할때 다른 항이 변하지 않으므로 변수로 따로 저장
			for (int j = k; j < _A.rows; j++) {
				_U.at[i][j] = _U.at[i][j] - (temp / _U.at[k][k]) * _U.at[k][j];
			}
		}
	}
}

void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P)
{
	copyVal(A, U);

	for (int k = 0; k < U.cols - 1; k++) {
		for (int i = k  + 1; i < U.rows ; i++) {
			double temp = U.at[i][k] / U.at[k][k];
			L.at[i][k] = temp;
			U.at[i][k] = 0;
			for (int j = k + 1; j < U.cols; j++) {
				U.at[i][j] = U.at[i][j] - temp * U.at[k][j];
			}
		}
	}
	for (int k = 0; k < A.cols; k++) {
		L.at[k][k] = 1;// L행렬에 identity matrix 더하기
	}
}

void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x)
{
	Matrix y = zeros(b.rows, 1);
	fwdsub(L, y, b);
	backSub(U, y, x);
}

void inv(Matrix A, Matrix Ainv)
{
	if (A.rows != A.cols) {
		printf("dimensin error");
		return;
	}
	Matrix _L = zeros(A.rows, A.cols);// nxn 행렬 생성
	Matrix _U = zeros(A.rows, A.cols);// nxn 행렬 생성
	Matrix _P = eye(A.rows);// identity 생성
	LUdecomp(A, _L, _U, _P);

	Matrix _b = zeros(A.rows, 1); //n x 1 행렬 생성
	Matrix _x = zeros(A.rows, 1); // n x 1 행렬 생성

	for (int i = 0; i <= Ainv.cols - 1; i++) {
		_b.at[i][0] = 1; //
		//if(i>0)
			//_b.at[i-1][0] = 0; // 초기화
		solveLU(_L, _U, _P, _b, _x);
		for (int j = 0; j <= Ainv.cols - 1; j++)
			Ainv.at[j][i] = _x.at[j][0];
		initMat(_b, 0);
		initMat(_x, 0);
	}
}


double func(double x) // Move this function to myNM.c
{
	/*double L = 4;
	double E = 70;
	double I = 52.9E-6;
	double w0 = 20;

	double F = w0 / (360 * L * E * I) * (7 * L * L * L * L - 30 * L * L * _x * _x + 15 * _x * _x * _x * _x);
	*/
	double F = 8 - 4.5 * (x - sin(x));
	return F;
}

double dfunc(double x)
{
	/*double L = 4;
	double E = 70;
	double I = 52.9E-6;
	double w0 = 20;

	double F = w0 / (360 * L * E * I) * (- 60 * L * L  * _x + 60  * _x * _x * _x);
	*/
	double F = -4.5 * (1 - cos(x));
	return F;
}



double bisectionNL(double _a0, double _b0, double _tol) // Move this function to myNM.c
{
	int k = 0;
	int Nmax = 1000;
	double a = _a0;
	double b = _b0;
	double xn = 0;
	double ep = 1000;

	do {
		xn = (a + b) / 2;
		ep = fabs(dfunc(xn));
		printf("Iteration:%d \t", k);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);

		if (dfunc(a) * dfunc(xn) < 0)
			b = xn;
		else
			a = xn;
		k++;
	} while (k<Nmax && ep>_tol);


	return xn;
}


Matrix  QRdecomp(Matrix _R) {
	int n = _R.rows;
	double normc = 0;
	Matrix I = eye(n);
	Matrix C = zeros(n, 1);
	Matrix E = zeros(n, 1);

	for (int count = 1; count < 44; count++) {
		Matrix Q = eye(_R.rows);
		for (int j = 0; j < _R.cols - 1; j++) {
			for (int i = 0; i < C.rows; i++) {
				C.at[i][0] = _R.at[i][j];
			}
			for (int k = 0; k <= j - 1; k++) {
				C.at[k][0] = 0;
			}
			normc = norm2(C);

			if (C.at[j][0] > 0)
				E.at[j][0] = normc;
			else
				E.at[j][0] = -1 * normc;

			Matrix V = addMat(C, E);
			Matrix V_tr = transpose(V);
			Matrix V_array = multipleMat(V, V_tr);
			Matrix V_value = multipleMat(V_tr, V);
			double d = 2 / V_value.at[0][0];
			Matrix V_d = smultipleMat(d, V_array);

			Matrix H = subMat(I, V_d);

			Q = multipleMat(Q, H);
			_R = multipleMat(H, _R);

			initMat(C, 0); initMat(E, 0); initMat(V_array, 0); initMat(H, 0);
			initMat(E, 0); initMat(V_tr, 0); initMat(V_value, 0); initMat(V_d, 0);
		}
		_R = multipleMat(_R, Q);
	}
	printMat(_R, "U");
	return _R;
	
}

Matrix eig(Matrix A)
{
	if (A.rows != A.cols) {
		printf("dimensin error");// Check Square Matrix
	}

	Matrix lamda = zeros(A.rows,1);
	Matrix U = QRdecomp(A);
		
	for (int m = 0; m<lamda.rows ; m++)
		lamda.at[m][0] = U.at[m][m];// Eigen Values
	return lamda;

	printMat(lamda, "Lamda");
		
}

double cond(Matrix A) {
	double cond = 0;
	double max = 0;
	double min = 0;
	
	Matrix Cond = eig(A);
	for (int i = 0; i < Cond.rows; i++) {
		if (sqrt(fabs(Cond.at[i][0])) > max) {
			max = sqrt(fabs(Cond.at[i][0]));
			if (min == 0)
				min = sqrt(fabs(Cond.at[i][0]));// min 값이 비워져 있을 경우 max값을 min값에 동일하게 입력
		}
		else if (sqrt(fabs(Cond.at[i][0])) < min) {
			if (Cond.at[i][0] < 0.000001)// 0으로 나눠지는 경우 배제
				min = min;
			else if (Cond.at[i][0] > 0.000001)
				min = sqrt(fabs(Cond.at[i][0]));
		}
		else
			;
	}
	cond = max / min;
	printf("Max.Eigenvalue is: %f  \n", max);
	printf("Max.Eigenvalue is: %f  \n", min);
	printf("Cond is: %f  \n", cond);
	return cond;
}

Matrix	linearFit(Matrix _x, Matrix _y) {
	int mx = _x.rows;
	int my = _y.rows;

	double a1 = 0;
	double a0 = 0;
	double DEN = 0;

	if ((mx != my) || (mx == 1))
		printf("ERROR: length of x and y must be equal and more than 1");
	else {
		double Sx = 0;
		double Sxx = 0;
		double Sxy = 0;
		double Sy = 0;

		for (int k = 0; k < mx; k++) {
			Sxx += _x.at[k][0] * _x.at[k][0];
			Sx += _x.at[k][0];
			Sy += _y.at[k][0];
			Sxy += _x.at[k][0] * _y.at[k][0];
		}

		DEN = mx * Sxx - Sx * Sx;
		a1 = (mx * Sxy - Sx * Sy) / DEN;
		a0 = (Sxx * Sy - Sxy * Sx) / DEN;

	}

	double z_array[] = { a1, a0 };
	return arr2Mat(z_array, 2, 1);
}

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}

void linearInterp(Matrix _x, Matrix _y, Matrix xq) {
	int mx = _x.rows;
	int my = _y.rows;

	for (int i = 0; i < xq.rows; i++) {
		int T = i * 5;

		if (i == 0) {
			xq.at[i][0] = _y.at[i][0];
			printf("When T is 0 C. Pressure is: %f\n", _y.at[i][0]);
		}

		else if (i % 2 == 1) {
			int r_i = ceil(i / 2); // T간격이 5일때 index를 표현하기 위함
			xq.at[i][0] = _y.at[r_i][0] * (T - _x.at[r_i + 1][0]) / (_x.at[r_i][0] - _x.at[r_i + 1][0]) + _y.at[r_i + 1][0] * (T - _x.at[r_i][0]) / (_x.at[r_i + 1][0] - _x.at[r_i][0]);
			printf("When T is %d C. Pressure is: %f\n", i * 5, xq.at[i][0]);
		}
		else {
			int index = i / 2;
			xq.at[i][0] = _y.at[index][0];
			printf("When T is %d C. Pressure is: %f\n", i * 5, xq.at[i][0]);
		}
	}
}

// Return the dy/dx results for the input data. (truncation error: O(h^2))
Matrix	gradient(Matrix _x, Matrix _y) {
	int n = _x.rows;
	int ny = _y.rows;

	//Check if n==ny
	if ((n != ny) || (n == 1))
		printf("ERROR: length of x and y must be equal and more than 1");
	else {
		Matrix Out = createMat(_x.rows, 1);

		//Assuming constant h
		double h = _x.at[1][0] - _x.at[0][0];

		//Assumption n>2
		if (n > 2) {
			// 1. 3-point FW
			Out.at[0][0] = (-3 * _y.at[0][0] + 4 * _y.at[1][0] - _y.at[2][0]) / (2 * h);

			// 2. 2-point central diff
			for (int i = 1; i < n - 1; i++) {
				Out.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);
			}

			// 3. 3-point BW
			Out.at[n - 1][0] = (_y.at[n - 3][0] - 4 * _y.at[n - 2][0] + 3 * _y.at[n - 1][0]) / (2 * h);

		}
		else
		{
			Out.at[0][0] = 0;//use 2-point forward
			Out.at[1][0] = 1;//use 2-point backward
		}
		return Out;
	}

}

void gradient1D(double x[], double y[], double dydx[], int m) {
	Matrix X = arr2Mat(x, m, 1);
	Matrix Y = arr2Mat(y, m, 1);
	Matrix Dydx = arr2Mat(dydx, m, 1);

	if ((X.rows != Y.rows) || (m == 1))
		printf("ERROR: length of x and y must be equal and more than 1");
	else {

		//Assuming constant h
		double h = X.at[1][0] - X.at[0][0];

		//Assumption n>2
		if (m > 2) {
			// 1. 3-point FW
			Dydx.at[0][0] = (-3 * Y.at[0][0] + 4 * Y.at[1][0] - Y.at[2][0]) / (2 * h);

			// 2. 2-point central diff
			for (int i = 1; i < m - 1; i++) {
				Dydx.at[i][0] = (Y.at[i + 1][0] - Y.at[i - 1][0]) / (2 * h);
			}

			// 3. 3-point BW
			Dydx.at[m - 1][0] = (Y.at[m - 3][0] - 4 * Y.at[m - 2][0] + 3 * Y.at[m - 1][0]) / (2 * h);

		}
		else
		{
			Dydx.at[0][0] = 0;//use 2-point forward
			Dydx.at[1][0] = 1;//use 2-point backward
		}
		printMat(Dydx, "Result:");
	}
}

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
Matrix	gradientFunc(double myFunc(const double x), Matrix xin) {
	int n = xin.rows;
	Matrix y = createMat(n, 1);

	//define y[0] to y[n-1]
	for (int i = 0; i < n; i++)
		y.at[i][0] = myFunc(xin.at[i][0]);

	//Use gradient() Numerical differentiation
	return gradient(xin, y);
}

double newtonRaphson(double _x0, double _tol)
{
	double xn = _x0;
	double xn1 = 0;
	double ep = 1000;
	int Nmax = 1000;
	int k = 0;

	do {
		xn1 = xn - func(xn) / dfunc(xn);
		ep = fabs((xn1 - xn) / xn);
		printf("Iteration:%d \t", k);
		printf("X(k): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);
		k++;
		xn = xn1;
	} while (k < Nmax && ep > _tol);

	return xn;
	printf("\n");
}

double newtonRaphsonFunc(double myFunc(const double x), double mydFunc(const double x), double x0, double tol)
{
	double xn = x0;
	double xn1 = 0;
	double ep = 1000;
	int Nmax = 1000;
	int k = 0;

	do {
		xn1 = xn - myFunc(xn) / mydFunc(xn);
		ep = fabs((xn1 - xn) / xn);
		printf("Iteration:%d \t", k);
		printf("X(k): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);
		k++;
		xn = xn1;
	} while (k < Nmax && ep > tol);

	return xn;
}

double trapz(double x[], double y[], int m) {
	int N = m - 1;
	double I = 0;
	double interval = x[1] - x[0];
	for (int i = 0; i < N; i++)
		I += (y[i+1]+ y[i]) * interval;

	I = I / 2;
	return I;
}

double integral(double func(const double x), double a, double b, int n) {
	double h = (b - a) / n;
	double f_x0 = func(a);
	double f_x1 = 0;
	double f_x2 = 0;
	double f_x3 = func(b);
	double I = 0;

	for (int i = 1; i <= n - 1; i= i + 2)
		f_x1 += func(a + i * h);
	f_x1 = f_x1 * 4;

	for (int k = 2; k <= n - 2; k= k + 2)
		f_x2 += func(a + k * h);
	f_x2 = f_x2 * 2;

	I = (f_x0+ f_x1+ f_x2+ f_x3)*h/3;
	return I;
}
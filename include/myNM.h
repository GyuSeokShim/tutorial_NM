/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include "myMatrix.h"

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);

// Matrix subtraction
extern	Matrix	subMat(Matrix _A, Matrix _B);

// Apply fwd-substitution
extern	void	fwdsub(Matrix _L, Matrix _y, Matrix _b);

// Apply back-substitution
extern	void	backSub(Matrix _U, Matrix _y, Matrix _x);

// Gauss-Elimination
extern	void	gaussElim(Matrix _A, Matrix _b, Matrix _U, Matrix _d);

// LU decomposition
extern	void	LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P);

// Function that solves for Ax=LUx=b, permutation P is applied
extern	void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x);
//extern	double	solveLU(Matrix L, Matrix U, Matrix P, Matrix b);

extern void inv(Matrix A, Matrix Ainv);

extern float bisectionNL(float _a0, float _b0, float _tol);

extern float func(float _x);

extern float dfunc(float _x);

extern float newtonRaphson(float _x0, float _tol);

extern	Matrix QRdecomp(Matrix _R);

extern	Matrix  eig(Matrix A);

extern double cond(Matrix A);

extern Matrix linearFit(Matrix _x, Matrix _y);

extern Matrix arr2Mat(double* _1Darray, int _rows, int _cols);

extern void linearInterp(Matrix _x, Matrix _y, Matrix xq);

extern Matrix	gradient(Matrix _x, Matrix _y);

extern void	gradient1D(double x[], double y[], double dydx[], int m);

extern Matrix	gradientFunc(double myFunc(const double x), Matrix xin);

extern double newtonRaphsonFunc(double myFunc(const double x), double mydFunc(const double x), double x0, double tol);

extern double trapz(double x[], double y[], int m);

extern double integral(double func(const double x), double a, double b, int n);

#endif
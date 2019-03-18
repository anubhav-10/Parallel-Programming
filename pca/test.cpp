#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#define epsilon 1e10

using namespace std;


void transpose(float *Data, int M, int N, float *Data_T) {
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			Data_T[j * M + i] = Data[i * N + j];
		}
	}
}
void matmul(float *mat1, float *mat2, int M, int N, int K, float *res) {
	for(int i = 0; i < M; i++){
		for(int j = 0; j < K; j++){
			for(int k = 0; k < N; k++)
				res[i * K + j] += mat1[i * N + k] * mat2[k * K + j];
		}
	}
}

float innerProduct(float *A, float *B, int N) {
	float iProd = 0;
	for(int i = 0; i < N; i++) {
		iProd += A[i] * B[i];
	}
	return iProd;
}

void projection(float *u, float *a, int N, float *res) {
	float scalar = innerProduct(u, a, N) / innerProduct(u, u, N);
	for(int i = 0; i < N; i++) {
		res[i] = u[i] * scalar;
	}
}

void copyVector(float *from, float *to, int N) {
	for(int i = 0; i < N; i++)
		to[i] = from[i];
}

void classicalGS(float *mat, int M, int N, float *Q, float *R) {
	// we need a_i's as columns so doing transpose
	// here M = N so a square matrix
	float mat_T[M * N];
	transpose(mat, M, N, mat_T);
	float E[M * N];
	float U[M * N];
	for(int i = 0; i < N; i++) {
		copyVector(&mat_T[i * N], &U[i * N], N);
		for(int j = 0; j < i; j++) {
			float proj[N];
			projection(&U[j * N], &mat_T[i * N], N, proj);
			for(int k = 0; k < N; k++)
				U[i * N + k] -= proj[k];
		}
		int mod = sqrt(innerProduct(&U[i * N], &U[i * N], N));
		for(int j = 0; j < N; j++) {
			E[i * N + j] = U[i * N + j] / mod;
		}
	}
	// Columns are still represented as rows
	transpose(E, N, N, Q);
	// matmul(E, mat, N, N, N, R);		
	// for(int i = 0; i < N; i++){
	// 	for(int j = 0; j < N; j++){
	// 		cout << Q[i * N + j] << " ";
	// 	}
	// 	cout<<endl;
	// }
	for(int i = 0; i < N; i++) {
		for(int j = i; j < N; j++) {
			R[i * N + j] = innerProduct(&E[i * N], &mat_T[j * N], N);
		}
	}

	// for(int i = 0; i < N; i++){
	// 	for(int j = 0; j < N; j++){
	// 		cout << R[i * N + j] << " ";
	// 	}
	// 	cout<<endl;
	// }
}

bool checkConvergence(float *D, float *E, float *D1, float *E1, int N) {
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(abs(D1[i * N + j] - D[i * N + j]) > epsilon || abs(E1[i * N + j] - E[i * N + j]) > epsilon)
				return 0;
		}
	}
	return 1;
}

void QRAlgo(float *mat, int M, int N, float *eigenValues, float *eigenVectors) {
	// float D[M * N];
	// float E[M * N] = {0};
	float *D = (float*)calloc(M * N, sizeof(float));
	float *E = (float*)calloc(M * N, sizeof(float));
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			D[i * N + j] = mat[i * N + j];
		}
		E[i * N + i] = 1;
	}
	bool converged = 0;
	while(!converged) {
		float Q[M * N], R[M * N] = {0};
		classicalGS(D, M, N, Q, R);
		float D1[M * N], E1[M * N];
		matmul(R, Q, N, N, N, D1);
		matmul(E, Q, N, N, N, E1);
		converged = checkConvergence(D, E, D1, E1, N);
		D = D1;
		E = E1;
	}
	for(int i = 0; i < N; i++) {
		eigenValues[i] = D[i * N + i];	
	}
	// it is colmuns of E
	eigenVectors = E;
	// free(D);
	// free(E);
}

int main(){
	float a[6] = {1,2,3,4,5,6};
	float b[4] = {1,2,3,4};
	float t[9] = {12, -51, 4, 6, 167, -68, -4, 24, -41};
	float Q[9], R[9];
	float M[4] = {24, -15, -15, 25};
	float eigenValues[2], eigenVectors[4];
	QRAlgo(M, 2, 2, eigenValues, eigenVectors);
	// classicalGS(t, 3, 3, Q, R);
	// matmul(a,b, 2,2,2, t);
	// transpose(a, 2, 3, Q);
	printf("%f %f\n", eigenValues[0], eigenValues[1]);
	// printf("%f %f\n", Q[2], Q[3]);
	// printf("%f %f\n", Q[4], Q[5]);
}
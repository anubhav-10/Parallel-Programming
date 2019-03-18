#include <stdio.h>
#include <iostream>
#include <cmath>
using namespace std;


void transpose(float *Data, int M, int N, float *Data_T){
	for(int i = 0; i < M; i++){
		for(int j = 0; j < N; j++){
			Data_T[j*N+i] = Data[i*M+j];
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

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			cout << R[i * N + j] << " ";
		}
		cout<<endl;
	}
}

int main(){
	float a[4] = {1,2,3,4};
	float b[4] = {1,2,3,4};
	float t[9] = {12, -51, 4, 6, 167, -68, -4, 24, -41};
	float Q[9], R[9];
	classicalGS(t, 3, 3, Q, R);
	// matmul(a,b, 2,2,2, t);
	// printf("%f %f\n", t[0], t[1]);
	// printf("%f %f\n", t[2], t[3]);
}

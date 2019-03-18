#include <malloc.h>
#include <omp.h>
#include <vector>
#include <cmath>
#define pb push_back
#define epsilon 1e5
using namespace std;

void transpose(float *Data, int M, int N, float *Data_T) {
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			Data_T[j*N+i] = Data[i*M+j];
		}
	}
}

void matmul(float *mat1, float *mat2, int M, int N, int K, float *res) {
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < K; j++) {
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
	for(int i = 0; i < N; i++) {
		for(int j = i; j < N; j++) {
			R[i * N + j] = innerProduct(&E[i * N], &mat_T[j * N], N);
		}
	}
}

bool checkConvergence(float *D, float *E, float *D1, float *E1, int N) {
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(abs(D1[i * N + j] - D[i * N + J]) > epsilon || abs(E1[i * N + j] - E[i * N + J]) > epsilon)
				return 0;
		}
	}
	return 1;
}

void QRAlgo(float *mat, int M, int N, float *eigenValues, float *eigenVectors) {
	float D[M * N];
	float E[M * N] = {0};
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
	eigenVectors = E;
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T) {
	float Data_T[M * N];
	transpose(D, M, N, Data_T);
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}

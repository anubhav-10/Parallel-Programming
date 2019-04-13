#include <malloc.h>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <algorithm>
#define pb push_back
#define epsilon 1e-3
using namespace std;

void print_matrix(string name, int M, int N, float* A){
	cerr << name << ": \n";
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			cerr << A[i*N + j] << " " ;
		}
		cerr << endl;
	}
	cerr << endl;
}

void transpose(float *Data, int M, int N, float *Data_T) {
	#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			Data_T[j * M + i] = Data[i * N + j];
		}
	}
}

void matmul(float *mat1, float *mat2, int M, int N, int K, float *res) {
	#pragma omp parallel for
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < K; j++) {
			res[j * M + i] = 0;
		}
	}

	#pragma omp parallel for collapse(2)
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < K; j++) {
			for(int k = 0; k < N; k++)
				res[i * K + j] += mat1[i * N + k] * mat2[k * K + j];
		}
	}
}

void inverse(float *mat, int M, int N, float *mat_inv) {
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			mat_inv[j * M + i] = 1 / mat[i * N + j];
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

void GS(float *mat, int M, int N, float *Q, float *R) {
	float mat_T[M * N];
	transpose(mat, M, N, mat_T);
	float V_T[M * N];
	float Q_T[M * N];
	// #pragma omp parallel for
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++)
			V_T[i * N + j] = mat_T[i * N + j];
		for(int j = 0; j < i; j++) {
			R[j * N + i] = innerProduct(&Q_T[j * N], &mat_T[i * N], N);
			for(int k = 0; k < N; k++) {
				V_T[i * N + k] -= R[j * N + i] * Q_T[j * N + k];
			}
		}
		R[i * N + i] = sqrt(innerProduct(&V_T[i * N], &V_T[i * N], N));
		for(int j = 0; j < N; j++) {
			Q_T[i * N + j] = V_T[i * N + j] / R[i * N + i];
		}
	}
	transpose(Q_T, N, N, Q);
}

bool checkConvergence(float *D, float *E, float *D1, float *E1, int N) {
	float conv = 0;
	// #pragma omp parallel for reduction(+:conv)
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			conv += abs(D1[i * N + j] - D[i * N + j]) + abs(E1[i * N + j] - E[i * N + j]);
			// if(abs(D1[i * N + j] - D[i * N + j]) > epsilon || abs(E1[i * N + j] - E[i * N + j]) > epsilon)
			// 	return 0;
		}
	}
	// cout << conv << endl;
	return conv < epsilon;
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
		GS(D, M, N, Q, R);
		float D1[M * N], E1[M * N];
		matmul(R, Q, N, N, N, D1);
		matmul(E, Q, N, N, N, E1);
		converged = checkConvergence(D, E, D1, E1, N);

		#pragma omp parallel for
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < N; j++) {
				D[i * N + j] = D1[i * N + j];
				E[i * N + j] = E1[i * N + j];
			}
		}

		// D = D1;
		// E = E1;
	}
	for(int i = 0; i < N; i++) {
		eigenValues[i] = D[i * N + i];	
	}
	// it is colmuns of E
	// print_matrix("eigenVectors", N, N, E);
	for(int i = 0; i < N * N; i++)
		eigenVectors[i] = E[i];
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T) {
	float D_T[N * M];
	transpose(D, M, N, D_T);
	float mat[N * N];
	matmul(D_T, D, N, M, N, mat);

	float eigenValues[N * N], eigenVectors[N * N];
	QRAlgo(mat, N, N, eigenValues, eigenVectors);

	// sort(eigenValues, eigenValues + N);
	float sigma[N * M];
	float sigma_inv[M * N] = {0};
	for(int i = 0; i < N; i++) {
		(*SIGMA)[i] = sqrt(eigenValues[i]);
		sigma[i * M + i] = sqrt(eigenValues[i]);
		sigma_inv[i * N + i] = 1 / sigma[i * M + i];
	}

 	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			(*U)[i * N + j] = eigenVectors[i * N + j];
		}
	}

	float U_T[N * N];
	transpose(*U, N, N, U_T);

	float temp[M * N];
	matmul(sigma_inv, U_T, M, N, N, temp);
	matmul(temp, D_T, M, N, M, *V_T);
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K) {
    float sum_sigma = 0;
    for(int i = 0; i < N; i++) {
            sum_sigma += SIGMA[i] * SIGMA[i];
    }
    float ret = (float)retention / 100.0;
    float variance = 0;
    *K = -1;
    // cout << sum_sigma << endl;
    // print_matrix("U", 1, N, U);
    for(int i = 0; i < N;i++) {
            variance += (SIGMA[i] * SIGMA[i]) / sum_sigma;
            if(variance > ret) {
                    *K = i + 1;
                    break;
            }
    }
    if(*K == -1) *K = N;

    *D_HAT = (float*)malloc(sizeof(float) * M * (*K));
    float W[N * (*K)];
	#pragma omp parallel for
    for(int i = 0; i < N; i++) {
            for(int j = 0; j < (*K); j++) {
                    W[i * (*K) + j] = U[i * N + j];
            }
    }
    matmul(D, W, M, N, *K, *D_HAT);
    // cerr << "K = " << *K << endl;
    // print_matrix("D-Hat", N, *K, *D_HAT);
}

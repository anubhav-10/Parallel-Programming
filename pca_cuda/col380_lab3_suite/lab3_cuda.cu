#include "lab3_cuda.h"
#include <iostream>
#include <fstream>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <algorithm>
using namespace std;

#define epsilon 1e-5
#define FILENAME "testcase_1000_300"

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void print_matrix(string name, int M, int N, double* A){
	cerr << name << ": \n";
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			cerr << A[i*N + j] << " " ;
		}
		cerr << endl;
	}
	cerr << endl;
}

void transpose(double *Data, int M, int N, double *Data_T) {
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < N; j++) {
			Data_T[j * M + i] = Data[i * N + j];
		}
	}
}

void matmul(double *mat1, double *mat2, int M, int N, int K, double *res) {
	for(int i = 0; i < M; i++) {
		for(int j = 0; j < K; j++) {
			res[j * M + i] = 0;
		}
	}

	for(int i = 0; i < M; i++) {
		for(int j = 0; j < K; j++) {
			for(int k = 0; k < N; k++)
				res[i * K + j] += mat1[i * N + k] * mat2[k * K + j];
		}
	}
}

__global__ void chessTournament(int N, int *p, int *q) {
	int tid = threadIdx.x;
	int i = blockIdx.x;
	int a, b;
	a = (tid + i) % (N - 1);
	if (tid != 0) {
		b = ((N - tid) + i - 1) % (N - 1);
	}
	else{
		b = N - 1;
	}

	p[i * (N / 2) + tid] = min(a, b);
	q[i * (N / 2) + tid] = max(a, b);	
}

__global__ void findSinCos(double *data, int* _p, int* _q, int N, double *sin, double *cos, int i) {
	int k = _p[i * (N/2) + blockIdx.x];
	int l = _q[i * (N/2) + blockIdx.x];
	double p = data[k * N + l];
	double y = (data[k * N + k] - data[l * N + l]) / 2.0;
	double d = fabs(y) + sqrt(p * p + y * y);
	double r = sqrt(p * p + d * d);
	double c, s;
	if(fabs(d) < 1e-2 && fabs(p) < 1e-2) {
		c = 1;
		s = 0;
	}
	else {
		c = d / r;
		s = p / r;
	}

	sin[blockIdx.x] = s;
	cos[blockIdx.x] = c;

	// printf("sin: %f cos: %f sq: %f \n",s, c, s*s + c*c);
}

__global__ void updateRow(double *data, int* _p, int* _q, int N, double *sin, double *cos, int i, double *out) {
	int tid = threadIdx.x;
	__shared__ double s;
	__shared__ double c;
	if (threadIdx.x == 0) {
		s = sin[blockIdx.x];
		c = cos[blockIdx.x];
	}
	__syncthreads();
	int p = _p[i * (N/2) + blockIdx.x];
	int q = _q[i * (N/2) + blockIdx.x];
	int id1 = p * N + tid;
	int id2 = q * N + tid;

	double ap = data[id1];
	double aq = data[id2];
	out[id1] = c * ap - s * aq; 
	out[id2] = s * ap + c * aq;
}

__global__ void updateCol(double *data, int* _p, int* _q, int N, double *sin, double *cos, int i, double *out) {
	int tid = threadIdx.x;
	__shared__ double s;
	__shared__ double c;
	if (threadIdx.x == 0) {
		s = sin[blockIdx.x];
		c = cos[blockIdx.x];
	}
	__syncthreads();
	int p = _p[i * (N/2) + blockIdx.x];
	int q = _q[i * (N/2) + blockIdx.x];
	int id1 = tid * N + p;
	int id2 = tid * N + q;
	
	double ap = data[id1];
	double aq = data[id2];
	out[id1] = c * ap - s * aq; 
	out[id2] = s * ap + c * aq;
}

bool convergence(double *data, double *new_data, int N) {
	double diff = 0;
	for(int i = 0; i < N * N; i++) {
		diff += fabs(data[i] - new_data[i]);
	}
	cout << "convergence: " << diff << endl;
	return diff < epsilon;
}

void jacobi(double *D, int N, double *EIGENVALUES, double *EIGENVECTOR) {
	double *data;
	double *old_data = (double*)malloc(sizeof(double) * N * N);
	double *new_data = (double*)malloc(sizeof(double) * N * N);
	// double *temp1;

	for(int i = 0; i < N * N; i++) {
		old_data[i] = D[i];
	}
	for(int i = 0; i < N; i++) {
		EIGENVECTOR[i * N + i] = 1;
	}
	double *EIGENVECTORCuda;
	cudaMalloc((void**)&data, sizeof(double) * N * N);
	cudaMalloc((void**)&EIGENVECTORCuda, sizeof(double) * N * N);
	cudaMemcpy(data, D, sizeof(double) * N * N, cudaMemcpyHostToDevice);
	cudaMemcpy(EIGENVECTORCuda, EIGENVECTOR, sizeof(double) * N * N, cudaMemcpyHostToDevice);
	int *p;
	int *q;
	cudaMalloc((void**)&p, sizeof(int) * (N - 1) * (N/2));
	cudaMalloc((void**)&q, sizeof(int) * (N - 1) * (N/2));

	chessTournament<<<N - 1, N/2>>>(N, p, q);
	cudaDeviceSynchronize();

	double *sin, *cos;
	cudaMalloc((void**)&sin, sizeof(double) * (N / 2));
	cudaMalloc((void**)&cos, sizeof(double) * (N / 2));

	bool converged = 0;

	while(!converged) {
		for(int i = 0; i < N - 1; i++) {
			findSinCos<<<N / 2, 1>>>(data, p, q, N, sin, cos, i);
			cudaDeviceSynchronize();
			updateRow<<<N / 2, N>>>(data, p, q, N, sin, cos, i, data);
			cudaDeviceSynchronize();
			updateCol<<<N / 2, N>>>(data, p, q, N, sin, cos, i, data);
			updateCol<<<N / 2, N>>>(EIGENVECTORCuda, p, q, N, sin, cos, i ,EIGENVECTORCuda);
			cudaDeviceSynchronize();
		}
		cudaMemcpy(new_data, data, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
		converged = convergence(old_data, new_data, N);
		for(int j = 0; j < N * N; j++)
			old_data[j] = new_data[j];
	}
	cudaMemcpy(EIGENVECTOR, EIGENVECTORCuda, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
	// print_matrix("EIGENVECTOR", 4, 4, EIGENVECTOR);
	// print_matrix("new_data", N, N, new_data);
	// cout<<N<<endl;
	for(int i = 0; i < N; i++) {
		EIGENVALUES[i] = fabs(new_data[i * N + i]);
		// cout<<EIGENVALUES[i]<<",";
	}
	sort(EIGENVALUES, EIGENVALUES + N);
	for(int i = 0; i < N; i++) {
		cout<<EIGENVALUES[i]<<",";	
	}

	cout<<endl;
	cudaFree(p);
	cudaFree(q);
	cudaFree(data);
	free(old_data);
	free(new_data);
}

void SVD_and_PCA (int M, 
        int N, 
        double* D, 
        double** U, 
        double** SIGMA, 
        double** V_T, 
        double** D_HAT, 
        int *K,
        int retention) {
    // SVD
    *U = (double*)malloc(sizeof(double) * N * N);
    *SIGMA = (double*)malloc(sizeof(double) * N);
    *V_T = (double*)malloc(sizeof(double) * M * M);

    double D_T[N * M];
    transpose(D, M, N, D_T);
    double DTD[N * N];
    matmul(D_T, D, N, M, N, DTD);

	// print_matrix("res", N, N, DTD);

    double EIGENVALUES[N];
    double EIGENVECTOR[N * N];
    memset(EIGENVECTOR, 0, sizeof(EIGENVECTOR[0]) * N * N);
    jacobi(DTD, N, EIGENVALUES, EIGENVECTOR);

 //    double sigma_inv[M * N];
 //    memset(sigma_inv, 0, sizeof(sigma_inv[0]) * M * N);
 //    for(int i = 0; i < N; i++) {
 //    	(*SIGMA)[i] = sqrt(EIGENVALUES[i]);
 //    	// sigma[i * M + i] = sqrt(EIGENVALUES[i]);
 //    	sigma_inv[i * N + i] = 1 / sqrt(EIGENVALUES[i]);
 //    }

 //    for(int i = 0; i < N; i++) {
 //    	for(int j = 0; j < N; j++) {
 //    		(*U)[i * N + j] = EIGENVECTOR[i * N + j];
 //    	}
 //    }

	// double U_T[N * N];
	// transpose(*U, N, N, U_T);

	// // double temp[M * N];
	// double *temp = (double*)malloc(sizeof(double) * M * N);
	// matmul(sigma_inv, U_T, M, N, N, temp);
	// matmul(temp, D_T, M, N, M, *V_T);

	// free(temp);

	// // PCA
 //    double sum_sigma = 0;
 //    for(int i = 0; i < N; i++) {
 //            sum_sigma += (*SIGMA)[i] * (*SIGMA)[i];
 //    }
 //    double ret = (double)retention / 100.0;
 //    double variance = 0;
 //    *K = -1;
 //    for(int i = 0; i < N;i++) {
 //            variance += ((*SIGMA)[i] * (*SIGMA)[i]) / sum_sigma;
 //            if(variance > ret) {
 //                    *K = i + 1;
 //                    break;
 //            }
 //    }
 //    if(*K == -1) *K = N;

 //    *D_HAT = (double*)malloc(sizeof(double) * M * (*K));
 //    double W[N * (*K)];
 //    for(int i = 0; i < N; i++) {
 //            for(int j = 0; j < (*K); j++) {
 //                    W[i * (*K) + j] = (*U)[i * N + j];
 //            }
 //    }
 //    matmul(D, W, M, N, *K, *D_HAT);
 //    cerr << "K = " << *K << endl;
 //    // print_matrix("D-Hat", N, *K, *D_HAT);    
}

void read_file(char* filename, int M, int N, double* A) {
    ifstream ifile;
    ifile.open(filename, ios::in);

    double tmp;
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++){
            ifile >> tmp;
            A[i * N + j] = tmp;
        }
    }

    ifile.close();
}

int main(){
	int N = 300;
	int M = 1000;
	double *D = (double*)malloc(sizeof(double) * (M) * N);
	// double *D_T = (double*)malloc(sizeof(double) * (N) * N);
	// double *res = (double*)malloc(sizeof(double) * (N) * N);
	read_file((char*)FILENAME, 1000, 300, D);
	// transpose(D, 4, 4, D_T);
	// matmul(D_T, D, 4, 4, 4, res);
	// print_matrix("res", 4, 4, res);
	// return 0;
	// double *EIGENVECTOR = (double*)malloc(sizeof(double) * N * N);
	// double * EIGENVALUES = (double*)malloc(sizeof(double) * N * N);
	// jacobi(res, N, EIGENVECTOR, EIGENVECTOR);
	double *U;
	double *SIGMA;
	double *V_T;
	double *D_HAT;
	int K;
	int retention = 80;

	SVD_and_PCA(M, N, D, &U, &SIGMA, &V_T, &D_HAT, &K, retention);
	// free(D);
	// free(SIGMA);
	// free(V_T);
	// free(D_HAT);
}
#include "lab3_cuda.h"
#include <iostream>
using namespace std;
// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */

__device__ void chessTournament(int N, int *p, int *q) {
	int tid = threadIdx.x;
	int i = blockIdx.x;
	p[i * N / 2 + tid] = (tid + i) % (N - 1);
	if (tid != 0) {
		q[i * N / 2 + tid] = ((N - tid) + i - 1) % (N - 1);
	}
	else{
		q[i * N / 2 + tid] = N - 1;
	}
}

void jacobi(double *D, int N, double *EIGENVALUES, double *EIGENVECTOR) {
	double *data;
	cudaMalloc((void**)&data, sizeof(double) * N * N);
	cudaMemcpy(data, D, sizeof(double) * N * N, cudaMemcpyHostToDevice);
	int *p;
	int *q;
	cudaMalloc((void**)&p, sizeof(int) * (N - 1) * N/2);
	cudaMalloc((void**)&q, sizeof(int) * (N - 1) * N/2);
	chessTournament<<<N - 1, N/2>>>(N, p, q);
	cudaDeviceSynchronize();

	int *p1 = (int*)malloc(sizeof(int) * (N - 1) * N/2);
	int *q1 = (int*)malloc(sizeof(int) * (N - 1) * N/2);
	cudaMemcpy(p1, p, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
	cudaMemcpy(p2, q, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
	for(int i = 0; i < N - 1; i++) {
		cout<<"["<<i<<"]"<<endl;
		for(int j = 0; j < N/2; j++)
			cout<<p1[i * (N - 1) + j]<< " " << q1[i * (N - 1) + j] << endl;
		cout<<endl;
	}	
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
    // write your code here
    double *EIGENVALUES = (double*)malloc(sizeof(double) * N * N);
    double *EIGENVECTOR = (double*)malloc(sizeof(double) * N * N);
}

int main(){
	double *D;
	int N = 8;
	double *EIGENVECTOR;
	double * EIGENVALUES;
	jacobi(D, N, EIGENVECTOR, EIGENVECTOR);
}
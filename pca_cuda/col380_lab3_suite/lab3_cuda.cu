#include "lab3_cuda.h"
#include <iostream>
#include <fstream>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <ctime>
#include <ratio>
#include <chrono>
using namespace std;

#define epsilon 1e-5
#define BLOCK_SIZE 16
// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
__global__ void Muld(double*, double*, int, int, double*);

void Mul(double* A, double* B, int hA, int wA, int wB, double* C) {
	int size; 
	double* Ad;
	size = hA * wA * sizeof(double);
	cudaMalloc((void**)&Ad, size);
	cudaMemcpy(Ad, A, size, cudaMemcpyHostToDevice);
	double* Bd;
	size = wA * wB * sizeof(double);
	cudaMalloc((void**)&Bd, size);
	cudaMemcpy(Bd, B, size, cudaMemcpyHostToDevice); 
	double* Cd;
	size = hA * wB * sizeof(double);
	cudaMalloc((void**)&Cd, size);
	dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
	dim3 dimGrid(wB / dimBlock.x, hA / dimBlock.y);
	Muld<<<dimGrid, dimBlock>>>(Ad, Bd, wA, wB, Cd);
	cudaMemcpy(C, Cd, size, cudaMemcpyDeviceToHost);
	cudaFree(Ad);
	cudaFree(Bd);
	cudaFree(Cd);
}

__global__ void Muld(double* A, double* B, int wA, int wB, double* C){
	int bx = blockIdx.x;
	int by = blockIdx.y;
	int tx = threadIdx.x;
	int ty = threadIdx.y;
	int aBegin = wA * BLOCK_SIZE * by;
	int aEnd   = aBegin + wA - 1;
	int aStep  = BLOCK_SIZE;
	int bBegin = BLOCK_SIZE * bx;
	int bStep  = BLOCK_SIZE * wB;
	double Csub = 0;
	for (int a = aBegin, b = bBegin;a <= aEnd;a += aStep, b += bStep) {
		__shared__ double As[BLOCK_SIZE][BLOCK_SIZE];
		__shared__ double Bs[BLOCK_SIZE][BLOCK_SIZE];
		As[ty][tx] = A[a + wA * ty + tx]; 
		Bs[ty][tx] = B[b + wB * ty + tx];
		__syncthreads();
		for (int k = 0; k < BLOCK_SIZE; ++k)
			Csub += As[ty][k] * Bs[k][tx];
		__syncthreads();
	}
	int c = wB * BLOCK_SIZE * by + BLOCK_SIZE * bx;
	C[c + wB * ty + tx] = Csub;
}

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
	double y = (data[l * N + l] - data[k * N + k]) / 2.0;
	double d = fabs(y) + sqrt(p * p + y * y);
	double r = sqrt(p * p + d * d);
	double c, s;
	if(fabs(d) < 1e-2 && fabs(p) < 1e-2) {
		c = 1;
		s = 0;
	}
	else {
		c = d / r;
		s = (fabs(y)/y)*(p / r);
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
	out[tid*N + p] = c * ap - s * aq; 
	out[tid*N + q] = s * ap + c * aq;
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
	int id1 = tid + N * p;
	int id2 = tid + N * q;
	
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

void jacobi(double *D_N, int N, double *EIGENVALUES, double *EIGENVECTOR) {

	double *D = D_N;
	int N1 = N;
	if(N % 2 == 1 ){
		D = (double*)calloc((N+1)*(N+1), sizeof(double));

		for(int i=0; i<N; i++){
			for(int j=0; j<N; j++){
				D[i*(N+1) + j] = D_N[i*N + j];
			}
		}
		N = N + 1;
		D[N*N - 1] = 1;
	}



	double *data;
	double *old_data = (double*)malloc(sizeof(double) * N * N);
	double *new_data = (double*)malloc(sizeof(double) * N * N);
	double *EIGENVECTOR_temp = (double*)calloc(N * N, sizeof(double));
	// double *temp1;

	for(int i = 0; i < N * N; i++) {
		old_data[i] = D[i];
	}
	for(int i = 0; i < N; i++) {
		EIGENVECTOR_temp[i * N + i] = 1;
	}
	double *EIGENVECTORCuda;
	cudaMalloc((void**)&data, sizeof(double) * N * N);
	cudaMalloc((void**)&EIGENVECTORCuda, sizeof(double) * N * N);
	cudaMemcpy(data, D, sizeof(double) * N * N, cudaMemcpyHostToDevice);
	cudaMemcpy(EIGENVECTORCuda, EIGENVECTOR_temp, sizeof(double) * N * N, cudaMemcpyHostToDevice);
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

	cudaMemcpy(EIGENVECTOR_temp, EIGENVECTORCuda, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
	// print_matrix("EIGENVECTOR", 4, 4, EIGENVECTOR);
	// print_matrix("new_data", N, N, new_data);
	// cout<<N<<endl;
	for(int i = 0; i < N1; i++) {
		EIGENVALUES[i] = fabs(new_data[i * N + i]);
		// cout<<EIGENVALUES[i]<<",";
	}

	double *EIGENVECTOR_T = (double*)calloc(N * N, sizeof(double));

	transpose(EIGENVECTOR_temp, N, N, EIGENVECTOR_T);

    vector<pair<double, int>> EVal;
    for(int i = 0; i < N1; i++) {
    	EVal.push_back(make_pair(EIGENVALUES[i], i));
    }

	sort(EIGENVALUES, EIGENVALUES + N);
	reverse(EIGENVALUES, EIGENVALUES + N);

    sort(EVal.begin(), EVal.end());
    reverse(EVal.begin(), EVal.end());
    for(int i = 0; i < N1; i++) {
    	int k = EVal[i].second;
    	for(int j = 0; j < N1; j++) {
    		EIGENVECTOR[j * N1 + i] = EIGENVECTOR_T[j * N + k];
    	}
    }

    N = N1;


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

    // double D_T[N * M];
    double D_T = (double*)malloc(sizeof(double) * M * N);
    transpose(D, M, N, D_T);
    // double DTD[N * N];
    double DTD = (double*)malloc(sizeof(double) * N * N);

    // double *D_cuda, *D_T_cuda, *DTD_cuda;
    // cudaMalloc((void**)&D_cuda, sizeof(double) * M * N);
    // cudaMalloc((void**)&D_T_cuda, sizeof(double) * N * M);
    // cudaMalloc((void**)&DTD_cuda, sizeof(double) * N * N);
    Mul(D_T, D, N, M, N, DTD);
    // cudaMemcpy(D, D_cuda, sizeof(double) * M * N, cudaMemcpyHostToDevice);
    // cudaMemcpy(D_T, D_T_cuda, sizeof(double) * N * M, cudaMemcpyHostToDevice);

    // Mul(D_T_cuda, D_cuda, N, M, N, DTD_cuda);
    // cudaMemcpy(DTD, DTD_cuda, sizeof(double) * N * N, cudaMemcpyDeviceToHost);

	// print_matrix("res", N, N, DTD);

    double *EIGENVALUES = (double*)malloc(sizeof(double) * N);
    double EIGENVECTOR = (double*)malloc(sizeof(double) * N * N);
    memset(EIGENVECTOR, 0, sizeof(EIGENVECTOR[0]) * N * N);
    jacobi(DTD, N, EIGENVALUES, EIGENVECTOR);

    double sigma_inv = (double*)malloc(sizeof(double) * M * N);;
    memset(sigma_inv, 0, sizeof(sigma_inv[0]) * M * N);
    for(int i = 0; i < N; i++) {
    	(*SIGMA)[i] = sqrt(EIGENVALUES[i]);
    	// sigma[i * M + i] = sqrt(EIGENVALUES[i]);
    	sigma_inv[i * N + i] = 1 / sqrt(EIGENVALUES[i]);
    }

    for(int i = 0; i < N; i++) {
    	for(int j = 0; j < N; j++) {
    		(*U)[i * N + j] = EIGENVECTOR[i * N + j];
    	}
    }

	double U_T = (double*)malloc(sizeof(double) * N * N);;
	transpose(*U, N, N, U_T);

	double *temp = (double*)malloc(sizeof(double) * M * N);
	Mul(sigma_inv, U_T, M, N, N, temp);
	Mul(temp, D_T, M, N, M, *V_T);

	// double *sigma_inv_cuda, *U_T_cuda, *temp_cuda, *V_T_cuda;
 //    cudaMalloc((void**)&sigma_inv_cuda, sizeof(double) * M * N);
 //    cudaMalloc((void**)&U_T_cuda, sizeof(double) * N * N);
 //    cudaMalloc((void**)&temp_cuda, sizeof(double) * M * N);
 //    cudaMalloc((void**)&V_T_cuda, sizeof(double) * M * M);

 //    cudaMemcpy(sigma_inv, sigma_inv_cuda, sizeof(double) * M * N, cudaMemcpyHostToDevice);
 //    cudaMemcpy(U_T, U_T_cuda, sizeof(double) * N * N, cudaMemcpyHostToDevice);

 //    Mul(sigma_inv_cuda, U_T_cuda, M, N, N, temp_cuda);
 //    Mul(temp_cuda, D_T_cuda, M, N, M, V_T_cuda);

 //    cudaMemcpy(*V_T, V_T_cuda, sizeof(double) * M * M, cudaMemcpyDeviceToHost);


	// free(temp);

	// PCA
    double sum_sigma = 0;
    for(int i = 0; i < N; i++) {
            sum_sigma += (*SIGMA)[i] * (*SIGMA)[i];
    }
    double ret = (double)retention / 100.0;
    double variance = 0;
    *K = -1;
    for(int i = 0; i < N;i++) {
            variance += ((*SIGMA)[i] * (*SIGMA)[i]) / sum_sigma;
            if(variance > ret) {
                    *K = i + 1;
                    break;
            }
    }
    if(*K == -1) *K = N;

    *D_HAT = (double*)malloc(sizeof(double) * M * (*K));
    double W[N * (*K)];
    for(int i = 0; i < N; i++) {
            for(int j = 0; j < (*K); j++) {
                    W[i * (*K) + j] = (*U)[i * N + j];
            }
    }

    print_matrix("sigma", 1, N, *SIGMA);
    print_matrix("U", N, N, *U);
    print_matrix("V_T", M, M, *V_T);

    Mul(D, W, M, N, *K, *D_HAT);
    // double *W_cuda, *D_HAT_cuda;
    // cudaMalloc((void**)&sigma_inv_cuda, sizeof(double) * M * N);
    // cudaMalloc((void**)&U_T_cuda, sizeof(double) * N * N);

    cerr << "K = " << *K << endl;
    print_matrix("D-Hat", M, *K, *D_HAT);    
}

void read_file(char* filename, int *num_samples, int *num_features, double** A) {
	std::chrono::high_resolution_clock::time_point t1, t2;
	t1 = std::chrono::high_resolution_clock::now();

    ifstream ifile;
    ifile.open(filename, ios::in);
    int M, N;
    ifile >> M >> N;
    cout << M << " " << N << endl;
    *A = (double *)malloc(sizeof(double)*M*N);
    num_samples[0] = M;
    num_features[0] = N;
    double tmp;
    for (int i=0; i<M; i++) {
        for (int j=0; j<N; j++){
            ifile >> tmp;
            *((*A) + i*N + j) = tmp;
        }
    }

    ifile.close();

	t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

  	std::cout << "File Reading Time " << time_span.count() << " seconds.\n";


}

int main(int argc, char **argv){
	// int N = 300;
	// int M = 1000;
	int M, N;
	char* filename = argv[1];
	double *D; //= (double*)malloc(sizeof(double) * (M) * N);
	// double *D_T = (double*)malloc(sizeof(double) * (N) * N);
	// double *res = (double*)malloc(sizeof(double) * (N) * N);
	read_file((char*)filename, &M, &N, &D);
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
	int retention = stoi(argv[2]);

	SVD_and_PCA(M, N, D, &U, &SIGMA, &V_T, &D_HAT, &K, retention);
	// free(D);
	// free(SIGMA);
	// free(V_T);
	// free(D_HAT);
}
#include <malloc.h>
#include <omp.h>
#include <vector>
#define pb push_back
using namespace std;

void computeTranspose(vector<vector<float>> Data, int M, int N){
	vector<vector<float>> Data_T;
	for(int i = 0; i < N; i++){
		vector<float> temp(M);
		Data_T.pb(temp);
	}

	for(int i = 0; i < M; i++){
		for(int j = 0; j < N; j++){
			Data_T[j][i] = Data[i][j];
		}
	}

	return Data_T;
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T) {
	vector<vector<float>> Data;
	for(int i = 0; i < M; i++){
		vector<float> temp;
		for(int j = 0; j < N; j++){
			temp.pb(D[i * M + j]);
		}
		Data.pb(temp);
	}

	vector<vector<float>> Data_T = computeTranspose(Data, M, N);

}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    
}

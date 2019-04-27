#include "lab4_mpi.h"

#include <malloc.h>
#include "mpi.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#define pb push_back 

using namespace std;
void print_matrix(string name, int M, int N, bool* A){
	cerr << name << ": \n";
	for(int i=0; i<M; i++){
		for(int j=0; j<N; j++){
			cerr << A[i*N + j] << " " ;
		}
		cerr << endl;
	}
	cerr << endl;
}
void witness(string s, int p, int *A) {
	A[0] = 0;
	for(int i = 1; i < p; i++) {
		for(int j = 0; j < 2 * p - 1; j++) {
			if (s[j] == s[i + j])
				continue;
			else if (s[j] != s[i + j]){
				A[i] = j;
				break;
			}
		}
	}
}


int duel(char* text, int n, string pattern, int *witness, int i, int j) {
	int k = witness[j - i];
	if (j + k >= n)
		return i;
	else if (text[j + k] != pattern[k])
		return i;
	return j;
}

vector<int> non_periodic_matching(char *text, int n, string pattern, int *witness) {
	int block_size = ceil(pattern.length() / 2.0);
	int num_blocks = n / block_size;

	int possible_matches[num_blocks];

	for(int i = 0; i < num_blocks; i++) {
		int start = i * block_size;
		int end = min((i + 1) * block_size, n);

		for(int j = start; j < end; j++) {
			start = duel(text, n, pattern, witness, start, j);
		}
		possible_matches[i] = start;
	}

	std::vector<int> matches;

	for(int i = 0; i < num_blocks; i++) {
		bool flag = 0;
		for(int j = 0; j < pattern.length(); j++) {
			if (possible_matches[i] + j >= n || pattern[j] != text[possible_matches[i] + j]){
				flag = 1;
				break;
			}
		}

		if(!flag){
			matches.pb(possible_matches[i]);
		}
	}

	return matches;
}
// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
vector<int> periodic_matching_util(char *text, int n, string pattern, int p) {

	string non_periodic_pattern = pattern.substr(0, 2 * p - 1);
	// cout<<non_periodic_pattern<<endl;
	int wit[p];

	witness(non_periodic_pattern, p, wit);

	// cerr<<"witness"<<endl;
	// for(int i = 0;i < p; i++)
	// 	cerr<<wit[i]<<" ";
	// cerr<<endl;

	vector<int> possible_matches = non_periodic_matching(text, n, non_periodic_pattern, wit);

	// for(auto u: possible_matches)
	// 	cout<<u<<" ";
	// cout<<endl;

	int l = pattern.length();
	string u = pattern.substr(0, p);
	string v = pattern.substr(0, l % p);

	string uuv = u + u + v;
	// cout<<uuv<<endl;

	bool *M = (bool*)calloc(n, sizeof(bool));
	for (int i = 0; i < possible_matches.size(); i++) {
		int start = possible_matches[i];
		if (start + uuv.length() > n) break;
		bool flag = 0;
		for(int j = 0; j < uuv.length(); j++) {
			if (text[start + j] != uuv[j]){
				flag = 1;
				break;
			}
		}

		if(!flag)
			M[possible_matches[i]] = 1;
	}

	// for(int i = 0; i < n; i++)
	// 	cerr<<M[i]<< " ";
	// cerr<<endl;

	int m = pattern.length();
	vector<int> matches;
	for(int i = 0; i < p; i++) {
		// int k = 0;
		int col = i;
		for(int j = i; j < n; j += p) {
			bool flag = 0;
			for(int k = 0; k < m/p - 1; k++){
				if (k * p + j >= n){
					flag = 1;
					break;
				}
				else if(!M[k * p + j]){
					flag = 1;
					break;
				}
			}
			if(!flag)
				matches.pb(j);
		}
	}

	free(M);

	return matches;
}

void periodic_pattern_matching (
		int n, 
		char *text, 
		int num_patterns, 
		int *m_set, 
		int *p_set, 
		char **pattern_set, 
		int **match_counts, 
		int **matches)
{	
	int num_proc;
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int start = rank * (num_patterns / num_proc);
	int end = (rank + 1) * (num_patterns / num_proc);
	if (rank == num_proc - 1)
		end = num_patterns;

	vector<int> matches_local;
	vector<int> matches_count_local;
	for(int i = start; i < end; i++) {
		vector<int> v;
		char *pat = pattern_set[i];
		string pattern(pat);
		v = periodic_matching_util(text, n, pattern, p_set[i]);
		matches_local.insert(matches_local.end(), v.begin(), v.end());
		matches_count_local.push_back(v.size());
	}

	//SB1:
	int *sbuf;
	if (matches_count_local.size()==0)
		sbuf = nullptr;
	else
		sbuf = &matches_count_local[0];
	int sendcnt = end - start;

	//RB1:
	// int *rb1;
	int *revcnts1 = (int*)malloc(sizeof(int) * num_proc);
	int *disp;
	if (rank == 0){
		for(int i = 0; i < num_proc; i++){
			revcnts1[i] = (num_patterns / num_proc);
		}
		if (num_proc > 1)
			revcnts1[num_proc - 1] = num_patterns - (num_proc - 2) * (num_patterns / num_proc);
		else {
			revcnts1[0] = num_patterns;
		}

		disp = (int*)malloc(sizeof(int) * num_proc);
		disp[0] = 0;
		for(int i = 1; i < num_proc; i++){
			disp[i] = disp[i - 1] + revcnts1[i - 1];
		}

		*match_counts = (int*)malloc(sizeof(int) * (disp[num_proc - 1] + revcnts1[num_proc - 1]));
	}

	MPI_Gatherv(sbuf, sendcnt, MPI_INT, *match_counts, revcnts1, disp, MPI_INT, 0, MPI_COMM_WORLD);

	// if(rank == 0){
	// 	for(int i=0;i<num_patterns;i++)
	// 		cout << (*match_counts)[i] <<" : ";
	// 	cout <<endl;
	// }
	int *sbuf2;
	if (matches_local.size() == 0)
		sbuf2 = nullptr;
	else
		sbuf2 = &matches_local[0];
	int sendcnt2 = matches_local.size();

	int *revcnts2 = (int*)calloc(num_proc, sizeof(int));
	int *disp2;
	int size;
	if (rank == 0) {
		for (int i = 0; i < num_proc; i++) {
			// revcnts2[i] += 
			int start = i * (num_patterns / num_proc);;
			int end = (i + 1) * (num_patterns / num_proc);

			if (i == num_proc - 1)
				end = num_patterns;
			revcnts2[i] = 0;
			// cout << start << " : "<< end << endl;

			for (int j = start; j < end; j++){
				revcnts2[i] += (*match_counts)[j];
			}
		}

		disp2 = (int*)malloc(sizeof(int) * num_proc);
		disp2[0] = 0;
		for(int i = 1; i < num_proc; i++){
			disp2[i] = disp2[i - 1] + revcnts2[i - 1];
		}

		size = disp2[num_proc - 1] + revcnts2[num_proc - 1];
		*matches = (int*)malloc(sizeof(int) * (disp2[num_proc - 1] + revcnts2[num_proc - 1]));
	}
	// if(rank == 0){
	// 	for(int i=0;i<num_proc;i++){
	// 		cout << revcnts2[i] << " : ";
	// 	}
	// 	cout << endl;

	// 	for(int i=0;i<num_proc;i++){
	// 		cout << disp2[i] << " : ";
	// 	}
	// 	cout << endl;
	

	// }

	MPI_Gatherv(sbuf2, sendcnt2, MPI_INT, *matches, revcnts2, disp2, MPI_INT, 0, MPI_COMM_WORLD);

	// if (rank == 0){
		
	// 	cout<<size<<endl;
	// 	sort(*matches, (*matches)+size);
	// 	for(int i = 0; i < size; i++)
	// 		cout << (*matches)[i] <<" ";
	// 	cout<<endl;
	// }
}

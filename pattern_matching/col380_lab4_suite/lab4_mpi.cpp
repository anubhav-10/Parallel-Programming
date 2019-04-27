#include "lab4_mpi.h"

#include <malloc.h>
// #include "mpi.h"
#include <vector>
#include <iostream>
#include <cmath>
#define pb push_back 

using namespace std;

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

	// bool *C = (bool*)calloc(p * (n / p), sizeof(bool));
	int m = pattern.length();
	vector<int> matches;
	for(int i = 0; i < p; i++) {
		int k = 0;
		int col = 0;
		for(int j = i; j < n; j += p) {
			if (M[j] == 1) k++;
			else k = 0;
			if (k >= (pattern.length()/p) - 1){
				// C[i * (n/p) + j] = 1;
				// cout<<j<<" "<<((k-1)*p)<<endl;
				// matches.pb(j - ((p-1) * p));
				// matches.pb(j - (m/p - 2)*p);
				matches.pb(i + j/p);
				// cout<<"here"<<endl;
				// break;
			}
			col++;
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

}

int main() {
	// char* text = "abcabcabcaabcabca";
	// string pattern = "abcabca";

	char* text = "babababababaabab";
	string pattern = "abababa";

	vector<int> matches = periodic_matching_util(text, 16, pattern, 2);
	for(auto u: matches)
		cout<<u<<" ";
	cout<<endl;
}
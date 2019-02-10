#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <climits>
#include <pthread.h>
#include "lab1_pthread.h"
#define pb push_back 
using namespace std;

class Point{
	public:
		float x, y, z;
		int cluster_id;
		Point() {}
		Point(float x, float y, float z){
			this->x = x;
			this->y = y;
			this->z = z;
			cluster_id = -1;
		}
		Point operator+ (const Point & p) const{
			return Point(x + p.x, y + p.y, z + p.z);
		}
		void print(){
			cerr<<x<<" "<<y<<" "<<z<<" "<<cluster_id<<endl;
		}
};

class Cluster{
	public:
		vector<Point> points;
		Point sum;
		int size;
		Point centroid;
		Cluster() {}
};

vector<Point> data;
vector<Cluster> cluster;
vector<Point> allCentroids;
int NUM_POINTS, NUM_CLUSTERS, NUM_THREADS, NUM_ITER;
bool flag;
pthread_mutex_t lock;
pthread_barrier_t centroid_barrier;
pthread_barrier_t convergence_barrier;

float distance(Point p1, Point p2){
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));	
}

const float MAX_ERROR = 1e-6;
bool convergence(Point p1, Point p2){
	return distance(p1, p2) < MAX_ERROR;
}

void* kmeans_thread(void *tid){
	int *id = (int*) tid;
	int len = NUM_POINTS/NUM_THREADS;
	int start = (*id)*len;
	int stop = start + len;
	if(*id == NUM_THREADS - 1)
		stop = NUM_POINTS;

	vector<Cluster> local_cluster(NUM_CLUSTERS);

	for(int x=0;x<1000;x++){
		for(int i=0;i<NUM_CLUSTERS;i++){
			local_cluster[i].size = 0;
			local_cluster[i].sum = Point(0, 0, 0);
		}

		for(int j=start;j<stop;j++){
			float dist = INT_MAX;
			Point p1 = data[j];
			int ID = 0;
			for(int k=0;k<NUM_CLUSTERS;k++){
				float d = distance(p1, cluster[k].centroid);
				if(d < dist){
					ID = k;
					dist = d;
				}
			}
			data[j].cluster_id = ID;
			local_cluster[ID].sum = local_cluster[ID].sum + p1;
			local_cluster[ID].size++;
		}

		pthread_mutex_lock(&lock);
		for(int i=0;i<NUM_CLUSTERS;i++){
			cluster[i].sum = cluster[i].sum + local_cluster[i].sum;
			cluster[i].size += local_cluster[i].size;
		}
		pthread_mutex_unlock(&lock);

		pthread_barrier_wait(&centroid_barrier);
		if(*id == 0){
			flag = 1;
			for(int i=0;i<NUM_CLUSTERS;i++){
				int s = cluster[i].size;
				Point p = cluster[i].sum;
				if(s!=0){
					p.x /= s;
					p.y /= s;
					p.z /= s;
					flag &= convergence(cluster[i].centroid, p);
					cluster[i].centroid = p;
					allCentroids.pb(p);
				}
			}
			NUM_ITER++;
			for(int j=0;j<NUM_CLUSTERS;j++){
				cluster[j].size = 0;
				cluster[j].sum = Point(0, 0, 0);
			}
		
		}
		pthread_barrier_wait(&convergence_barrier);
		if(flag)
			break;
	}
	return NULL;
}

void kmeans_pthread(int num_threads, int N, int K, int *data_points, int** data_point_clusters, float **centroids, int* num_iterations) {
	NUM_POINTS = N;
	NUM_CLUSTERS = K;
	NUM_THREADS = num_threads;
	NUM_ITER = 0;

	for(int i=0;i<N;i++){
		Point temp(*(data_points + i*3), *(data_points + i*3 + 1), *(data_points + i*3 + 2));
		data.pb(temp);
	}

	for(int i=0;i<K;i++){
		Cluster c;
		c.centroid = data[i];
		cluster.pb(c);
		allCentroids.pb(c.centroid);
	}

	pthread_t thr[num_threads];
	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&centroid_barrier, NULL, num_threads);
	pthread_barrier_init(&convergence_barrier, NULL, num_threads);

	int *tid = (int*) malloc (sizeof(int) * num_threads);
	*num_iterations = 0;
		// clear points vector
	for(int j=0;j<K;j++){
		cluster[j].size = 0;
		cluster[j].sum = Point(0, 0, 0);
	}

	// assign points to cluster
	for(int j=0;j<num_threads;j++){
		tid[j] = j;
		pthread_create(&thr[j], NULL, kmeans_thread, &tid[j]);
	}
	for(int j=0;j<num_threads;j++){
		pthread_join(thr[j], NULL);
	}
	pthread_barrier_destroy(&centroid_barrier);
	pthread_barrier_destroy(&convergence_barrier);
	*num_iterations = NUM_ITER;

	cout<<*num_iterations<<endl;
	*data_point_clusters = (int*)malloc(sizeof(int) * N*4);
	*centroids = (float*)malloc(sizeof(float) * (*num_iterations + 1)*K*3);
	for(int i=0;i<N;i++){
		*(*data_point_clusters + i*4) = (int)data[i].x;
		*(*data_point_clusters + i*4 + 1) = (int)data[i].y;
		*(*data_point_clusters + i*4 + 2) = (int)data[i].z;
		*(*data_point_clusters + i*4 + 3) = (int)data[i].cluster_id;
	}
	for(int i=0;i<*num_iterations+1;i++){
		for(int j=0;j<K;j++){
			*(*centroids + (i * K + j) * 3) = allCentroids[i * K + j].x;
			*(*centroids + (i * K + j) * 3 + 1) = allCentroids[i * K + j].y;
			*(*centroids + (i * K + j) * 3 + 2) = allCentroids[i * K + j].z;
		}
	}
}
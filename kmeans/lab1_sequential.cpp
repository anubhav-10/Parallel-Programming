#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <climits>
#include "lab1_sequential.h"
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
		// vector<Point> points;
		Point sum;
		int size;
		Point centroid;
		Cluster() {}
};

float distance(Point p1, Point p2){
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));	
}

const float MAX_ERROR = 1e-6;
bool convergence(Point p1, Point p2){
	return distance(p1, p2) < MAX_ERROR;
}

void kmeans_sequential(int N, int K, int *data_points, int **data_point_clusters, float **centroids, int *num_iterations){
	vector<Point> data;
	vector<Cluster> cluster;
	vector<Point> allCentroids;

	for(int i=0;i<N;i++){
		Point temp(*(data_points + i*3), *(data_points + i*3 + 1), *(data_points + i*3 + 2));
		data.pb(temp); 
	}
	// initialization using first K points
	// srand(time(0));
	for(int i=0;i<K;i++){
		Cluster c;
		c.centroid = data[i];
		cluster.pb(c);
		allCentroids.pb(c.centroid);
	}
	// total iterations
	*num_iterations = 0;
	for(int i=0;i<1000;i++){
		// clear points vector
		for(int j=0;j<K;j++){
			cluster[j].size = 0;
			cluster[j].sum = Point(0, 0, 0);
			// cluster[j].points.clear();
		}

		// assign points to cluster
		for(int j=0;j<data.size();j++){
			float dist = INT_MAX;
			Point p1 = data[j];
			int id = 0;
			for(int k=0;k<K;k++){
				float d = distance(p1, cluster[k].centroid);
				if(d < dist){
					id = k;
					dist = d;
				}
			}
			data[j].cluster_id = id;
			cluster[id].sum = cluster[id].sum + p1;
			cluster[id].size++;
			// cluster[id].points.pb(p1);
		}

		// recompute centroid
		bool flag = 1;
		for(int j=0;j<K;j++){
			int s = cluster[j].size;
			Point p = cluster[j].sum;
			// Point p(0, 0, 0);
			// for(int k=0;k<s;k++){
			// 	p.x += cluster[j].points[k].x;
			// 	p.y += cluster[j].points[k].y;
			// 	p.z += cluster[j].points[k].z;
			// }
			// cout<<s<<endl;
			if(s!=0){
				p.x /= s;
				p.y /= s;
				p.z /= s;
				flag &= convergence(cluster[j].centroid, p);
				cluster[j].centroid = p;
				allCentroids.pb(p);
			}
		}
		*num_iterations+=1;
		if (flag)
			break;
	}
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

#ifndef DENOISING_H
#define DENOISING_H

#include "util.h"
#include "anisotropic_neighborhood.h"

enum FEATURE { FLAT, EDGE, CORNER };

void prepare_point_cloud_data(PointCloud & point_cloud, int k, PCDNumber r_r, vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods,
	vector<vector<pair<size_t, PCDNumber>>> & big_neighborhoods, vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods,
	vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals, PCDNumber & avg_distance, GEO::NearestNeighborSearch_var & NN);

vector<PointCloud::Normal> normal_filtering(PointCloud & mesh, int n_ns, PCDNumber sigma_nn, PCDNumber sigma_ns_ratio,
	vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals,
	PCDNumber & avg_distance, bool use_mesh_normals, vector<vector<pair<size_t, PCDNumber>>> & membership_functions);

FEATURE farthest_points(vector<PointCloud::Point> & points, bool use_avg, PCDNumber edge_threshold, vector<PointCloud::Point> & cluster_centroids,
	vector<int> & cluster_membership);

void normal_correction(PointCloud & point_cloud, vector<PointCloud::Normal> &normals, vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods,
	vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals,
	PCDNumber & avg_distance, GEO::NearestNeighborSearch_var & NN, PCDNumber tau_n, int n_nc, vector<vector<PointCloud::Normal>> & multi_normals);

void update_point_positions(PointCloud & point_cloud, PointCloud & original_point_cloud, vector<PointCloud::Normal> & normals,
	vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals,
	PCDNumber & avg_distance, GEO::NearestNeighborSearch_var & NN, int n_fp, PCDNumber upsilon_f,
	PCDNumber upsilon_ec, PCDNumber tau_o_ratio, vector<FEATURE> & classification, vector<PointCloud::Normal> & feature_dirs,
	vector<vector<PointCloud::Normal>> & multi_normals);

struct DenoisingParameters
{
	int n_ext = 9; //number of external iterations
	int n_int = 1; //number of internal iterations

	double alpha = 1.0; //alpha
	double beta = 0.1; //beta
	double gamma = 0.5; //gamma
	double a_0_ratio = 0.35;

	int k = 50; //number of nearest neighbors
	double r_r = 3.0; //regular neighborhood ratio
	double tau_n = 0.2; //normal difference threshold

	int n_ns = 7; //number of normal smoothing iterations
	PCDNumber sigma_nn = 0.3; //normal difference sigma for normal smoothing
	PCDNumber sigma_ns_ratio = 1.5; //spatial difference sigma ratio for normal smoothing (based on average edge length)
	bool use_mesh_normals = false;
	int n_nc = 2;

	double delta_cc_ratio = 0.5; //displacement ratio for more external feature line (based on average edge length)
	double theta_half = 0.959931; //theta/2: theta is the narrow neighborhood angle parameter ... 0.959931 = 55 degrees

	int n_fp = 3; // number of (local) iterations for flat point update
	double upsilon_f = 0.3; //step size for flat point update
	double upsilon_ec = 0.5; //step size for edge and corner points update
	double tau_o_ratio = 2.0; //ratio used to define the maximum point displacement from original positions (based on average edge length)

	void show()
	{
		cout << "**************************************************" << endl;
		cout << "PARAMETERS" << endl;
		cout << "**************************************************" << endl;
		cout << "n_ext: " << n_ext << endl;
		cout << "n_int: " << n_int << endl;

		cout << "alpha : " << alpha << endl;
		cout << "beta : " << beta << endl;
		cout << "gamma : " << gamma << endl;
		cout << "a_0_ratio : " << a_0_ratio << endl;

		cout << "k : " << k << endl;
		cout << "r_r : " << r_r << endl;
		cout << "tau_n : " << tau_n << endl;

		cout << "n_ns : " << n_ns << endl;
		cout << "sigma_nn : " << sigma_nn << endl;
		cout << "sigma_ns_ratio : " << sigma_ns_ratio << endl;
		cout << "n_nc : " << n_nc << endl;

		cout << "delta_cc_ratio : " << delta_cc_ratio << endl;
		cout << "theta_half : " << theta_half << endl;

		cout << "n_fp : " << n_fp << endl;
		cout << "upsilon_f : " << upsilon_f << endl;
		cout << "upsilon_ec : " << upsilon_ec << endl;
		cout << "tau_o_ratio : " << tau_o_ratio << endl;
		cout << "**************************************************" << endl;
	}
};

struct DenoisingData
{
	vector<vector<pair<size_t, PCDNumber>>> regular_neighborhoods;
	vector<vector<pair<size_t, PCDNumber>>> big_neighborhoods;
	vector<vector<pair<size_t, PCDNumber>>> small_neighborhoods;
	vector<PCDNumber> areas;
	vector<PointCloud::Normal> regular_normals;
	PCDNumber avg_distance;
	GEO::NearestNeighborSearch_var NN;

	vector<PointCloud::Normal> normals;

	vector<vector<pair<size_t, PCDNumber>>> adaptive_kernels;

	vector<pair<size_t, PCDNumber>> single_adaptive_kernel;

	vector<FEATURE> feature_classification;
	vector<PointCloud::Normal> feature_dirs;
	vector<vector<PointCloud::Normal>> multi_normals;

	
};

PointCloud denoise(PointCloud & point_cloud, DenoisingParameters & denoising_parameters, DenoisingData & denoising_data);

#endif // DENOISING_H
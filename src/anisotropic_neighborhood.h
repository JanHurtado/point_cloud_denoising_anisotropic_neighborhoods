#ifndef ANISOTROPIC_NEIGHBORHOOD_H
#define ANISOTROPIC_NEIGHBORHOOD_H


#include "util.h"
#include "solver.h"

void compute_anisotropic_neighborhood(PointCloud & point_cloud, int point_index, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals,
	vector<pair<size_t, PCDNumber>> & regular_neighborhood, vector<pair<size_t, PCDNumber>> & small_neighborhood, PCDNumber alpha, PCDNumber beta,
	PCDNumber gamma, PCDNumber a_0_ratio, vector<pair<size_t, PCDNumber>> & membership_function);

void compute_all_anisotropic_neighborhoods(PointCloud & point_cloud, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals,
	vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods, vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods,
	PCDNumber alpha, PCDNumber beta, PCDNumber gamma, PCDNumber a_0_ratio, vector<vector<pair<size_t, PCDNumber>>> & membership_functions);


#endif // ANISOTROPIC_NEIGHBORHOOD_H

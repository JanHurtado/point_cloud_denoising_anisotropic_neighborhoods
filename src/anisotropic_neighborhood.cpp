#include "anisotropic_neighborhood.h"

void compute_anisotropic_neighborhood(PointCloud & point_cloud, int point_index, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals, 
	vector<pair<size_t, PCDNumber>> & regular_neighborhood, vector<pair<size_t, PCDNumber>> & small_neighborhood, PCDNumber alpha, PCDNumber beta, 
	PCDNumber gamma, PCDNumber a_0_ratio, vector<pair<size_t, PCDNumber>> & membership_function)
{
	membership_function.clear();
	PCDNumber t_inf = MAX_PCDNUMBER;

	int num_points = regular_neighborhood.size();
	int num_points_small_neighborhood = small_neighborhood.size();
	PCDMatrix error_matrix(num_points, num_points);
	error_matrix.fill(0);
	PCDMatrix distance_matrix(num_points, num_points);
	distance_matrix.fill(0);
	PCDMatrix gradient_matrix(num_points, num_points);
	gradient_matrix.fill(0);
	PCDMatrix area_matrix(num_points, num_points);
	area_matrix.fill(0);
	PCDMatrix central_point_normal_difference_matrix(num_points, num_points);
	central_point_normal_difference_matrix.fill(0);
	PCDMatrix linear_term(num_points, 1);
	linear_term.fill(0);

	PCDNumber total_area = 0.0f;

	PointCloud::Point p = point_cloud.point(PointCloud::VertexHandle(point_index));

	for (size_t i = 0; i < num_points; i++)
	{
		size_t idx1 = regular_neighborhood[i].first;
		PCDNumber total_edge_length = 0.0f;
		for (size_t j = 0; j < num_points; j++)
		{
			size_t idx2 = regular_neighborhood[j].first;
			PCDNumber normal_dif = (regular_normals[idx1].normalize() - regular_normals[idx2].normalize()).length();
			error_matrix(i, j) = normal_dif;
			error_matrix(j, i) = normal_dif;
		}
		PCDNumber area = areas[idx1];
		PCDNumber distance = regular_neighborhood[i].second;
		PointCloud::Point neigh_p = point_cloud.point(PointCloud::VertexHandle(idx1));
		PointCloud::Point neigh_n = regular_normals[idx1];
		PCDNumber distance_to_plane = point_to_plane_distance(p,neigh_n,neigh_p);
		distance_matrix(i, i) = distance;
		area_matrix(i, i) = area;
		PCDNumber central_point_normal_difference = (regular_normals[point_index].normalize() - regular_normals[idx1].normalize()).length();
		central_point_normal_difference_matrix(i, i) = central_point_normal_difference;
		linear_term(i, 0) = (beta * ((distance+distance_to_plane) * areas[point_index] * areas[idx1])) + (gamma* (central_point_normal_difference * areas[point_index] * areas[idx1]));

		total_area += area;
	}

	PCDMatrix H(num_points, num_points);
	H.fill(0);
	PCDMatrix NH(num_points, num_points);
	NH.fill(0);
	PCDMatrix error_matrix_t = error_matrix.transpose();
	PCDMatrix gradient_matrix_t = gradient_matrix.transpose();

	H = alpha * (area_matrix * (error_matrix)* area_matrix) + gamma * (gradient_matrix_t * gradient_matrix);

	NH = (H + H.transpose())*0.5f;

	PCDNumber area_constraint = total_area * a_0_ratio;

	vector<PCDNumber> solution;
	quadratic_programming_solver(num_points, linear_term, NH, area_matrix, area_constraint, solution);
	for (int j = 0; j < solution.size(); j++)
	{
		membership_function.push_back(pair<size_t, PCDNumber>(regular_neighborhood[j].first, solution[j]));
	}
}


void compute_all_anisotropic_neighborhoods(PointCloud & point_cloud, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals, 
	vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods, vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods, 
	PCDNumber alpha, PCDNumber beta, PCDNumber gamma, PCDNumber a_0_ratio, vector<vector<pair<size_t, PCDNumber>>> & membership_functions)
{
	membership_functions.clear();
	membership_functions.resize(point_cloud.n_vertices());

#pragma omp parallel
#pragma omp for
	for (int v_idx = 0; v_idx < point_cloud.n_vertices(); v_idx++)
	{
		compute_anisotropic_neighborhood(point_cloud, v_idx, areas, regular_normals, regular_neighborhoods[v_idx], 
			small_neighborhoods[v_idx], alpha, beta, gamma, a_0_ratio, membership_functions[v_idx]);
	}
}
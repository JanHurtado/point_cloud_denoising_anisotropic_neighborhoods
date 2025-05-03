#include "denoising.h"

void prepare_point_cloud_data(PointCloud & point_cloud, int k, PCDNumber r_r, vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods,
	vector<vector<pair<size_t, PCDNumber>>> & big_neighborhoods, vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods,
	vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals, PCDNumber & avg_distance, GEO::NearestNeighborSearch_var & NN)
{
	regular_neighborhoods.clear();
	big_neighborhoods.clear();
	small_neighborhoods.clear();
	areas.clear();
	regular_normals.clear();
	int n_vert = point_cloud.n_vertices();
	regular_neighborhoods.resize(n_vert, vector<pair<size_t, PCDNumber>>());
	big_neighborhoods.resize(n_vert, vector<pair<size_t, PCDNumber>>());
	vector<vector<pair<size_t, PCDNumber>>> k_neighborhoods(n_vert, vector<pair<size_t, PCDNumber>>());
	vector<vector<pair<size_t, PCDNumber>>> k_filtered_neighborhoods(n_vert, vector<pair<size_t, PCDNumber>>());
	vector<vector<pair<size_t, PCDNumber>>> proj_neighborhoods(n_vert, vector<pair<size_t, PCDNumber>>());
	vector<vector<pair<size_t, PCDNumber>>> rings(n_vert, vector<pair<size_t, PCDNumber>>());
	small_neighborhoods.resize(n_vert, vector<pair<size_t, PCDNumber>>());
	areas.resize(n_vert);
	regular_normals.resize(n_vert);

	int min_n_neighbors = 20;

	avg_distance = 0;
	int edge_count = 0;

	PCDNumber epsilon_d = 0.01;

	PCDNumber r_s = 1.5;

	PCDNumber r_b = 2.5;

	/* 
	****************************************************************************************************************************************************
	Initializing knn structure 
	****************************************************************************************************************************************************
	*/

	GEO::initialize();
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("algo");
	GEO::CmdLine::import_arg_group("co3ne");

	double * vdata = new double[point_cloud.n_vertices() * 3];

	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		vdata[v_it->idx() * 3] = p[0];
		vdata[v_it->idx() * 3 + 1] = p[1];
		vdata[v_it->idx() * 3 + 2] = p[2];
	}

	NN = GEO::NearestNeighborSearch::create(3, "CNN");
	NN->set_points(point_cloud.n_vertices(), vdata);

	/*
	****************************************************************************************************************************************************
	Computing rough average distance, initial normals and knn-based neighborhoods
	****************************************************************************************************************************************************
	*/

	PCDNumber rough_avg_distance = 0;
	PCDNumber rough_avg_distance_count = 0;

	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		PCDNumber p_ptr[3]; p_ptr[0] = p[0]; p_ptr[1] = p[1]; p_ptr[2] = p[2];
		std::vector<GEO::index_t> neighs(k);
		std::vector<PCDNumber> sq_dists(k);
		NN->get_nearest_neighbors(k, p_ptr, neighs.data(), sq_dists.data());
		PointCloud::Point t_mean(0, 0, 0);
		PCDMatrix points(k, 3);
		for (int i = 0; i < neighs.size(); i++)
		{
			PointCloud::Point neigh_point = point_cloud.point(PointCloud::VertexHandle(neighs[i]));
			k_neighborhoods[v_it->idx()].push_back(pair<size_t, PCDNumber>(neighs[i], (p - neigh_point).norm()));
			t_mean += neigh_point;
			points(i, 0) = neigh_point[0];
			points(i, 1) = neigh_point[1];
			points(i, 2) = neigh_point[2];
			if (i > 0 && i < 7) // closest 6 points filtering used to compute rough average distance
			{
				rough_avg_distance += sqrt(sq_dists[i]);
				rough_avg_distance_count += 1.0;
			}
		}
		t_mean /= (double)k;
		for (int i = 0; i < neighs.size(); i++)
		{
			points(i, 0) = points(i, 0) - t_mean[0];
			points(i, 1) = points(i, 1) - t_mean[1];
			points(i, 2) = points(i, 2) - t_mean[2];
		}
		PCDMatrix cov = (points.transpose()*points) / (double)k;

		Eigen::JacobiSVD<PCDMatrix> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
		auto U = svd.matrixU();

		PointCloud::Normal normal(U(0, 2), U(1, 2), U(2, 2));
		normal = normal.normalize_cond();
		PointCloud::Normal prevNormal = point_cloud.normal(*v_it);
		prevNormal = prevNormal.normalize_cond();
		//Normal orientation consistency propagation
		PCDNumber dp = normal | prevNormal;
		if (dp < -0.1)
		{
			normal = -normal;
		}
		regular_normals[v_it->idx()] = normal;
	}
	rough_avg_distance /= rough_avg_distance_count;

	/*
	****************************************************************************************************************************************************
	Filtering knn neighbors to avoid the inclusion of points with normal following an opposite direction to the evaluated point. E.g. in a thin metal
	sheet surface knn neighborhoods can include points of the corresponding opposite solid face. We filter these points considering the normal 
	initialization. If an isolated point was found, e.g. a single point with wrong normal direction, we flip its normal. 
	This is implementation is considered for practical purposes.
	****************************************************************************************************************************************************
	*/
	for (size_t i = 0; i < k_neighborhoods.size(); i++)
	{
		PointCloud::Normal ni = regular_normals[i];
		for (size_t j = 0; j < k_neighborhoods[i].size(); j++)
		{
			PointCloud::Normal nj = regular_normals[k_neighborhoods[i][j].first];
			if ((ni | nj) > -0.1)
			{
				k_filtered_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_neighborhoods[i][j].first, k_neighborhoods[i][j].second));
			}
		}
		if (k_filtered_neighborhoods[i].size() < 3)
		{
			k_filtered_neighborhoods[i].clear();
			regular_normals[i] = -regular_normals[i];
			ni = regular_normals[i];
			for (size_t j = 0; j < k_neighborhoods[i].size(); j++)
			{
				PointCloud::Normal nj = regular_normals[k_neighborhoods[i][j].first];
				if ((ni | nj) > -0.1)
				{
					k_filtered_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_neighborhoods[i][j].first, k_neighborhoods[i][j].second));
				}
			}
		}
	}

	/*
	****************************************************************************************************************************************************
	Filtering points to be considered for Delaunay triangulation. Avoids the inclusion of points that are very close to the evaluated point in order to
	create a more reliable triangulation for the following operations.
	****************************************************************************************************************************************************
	*/
	PCDNumber rough_min_distance = rough_avg_distance * epsilon_d;
	for (size_t i = 0; i < k_filtered_neighborhoods.size(); i++)
	{
		proj_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][0].first, k_filtered_neighborhoods[i][0].second));
		for (size_t j = 1; j < k_filtered_neighborhoods[i].size(); j++)
		{
			if (k_filtered_neighborhoods[i][j].second > rough_min_distance)
			{
				proj_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][j].first, k_filtered_neighborhoods[i][j].second));
			}
		}
	}

	/*
	****************************************************************************************************************************************************
	Computing "one ring" neighborhoods (N_1) using local 2D Delaunay triangulations, global average distance, maximum distances in N_1's, and
	per-point areas.
	****************************************************************************************************************************************************
	*/
	vector<PCDNumber> max_distances_in_rings(point_cloud.n_vertices(), 0);
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		size_t idx = v_it->idx();
		PointCloud::Normal normal = regular_normals[idx];
		size_t ns = proj_neighborhoods[idx].size();
		double * points2D = new double[ns * 2];
		PointCloud::Point rnp;
		point_to_plane_distance(point_cloud.point(PointCloud::VertexHandle(proj_neighborhoods[idx][2].first)), normal, p, rnp);
		for (int i = 0; i < proj_neighborhoods[idx].size(); i++)
		{
			PointCloud::Point neigh_point = point_cloud.point(PointCloud::VertexHandle(proj_neighborhoods[idx][i].first));
			PointCloud::Point vec = neigh_point - p;
			PointCloud::Point proj_neigh_point;
			point_to_plane_distance(neigh_point, normal, p, proj_neigh_point);
			PointCloud::Point centered_proj_neigh_point = proj_neigh_point - p;
			PointCloud::Point u = (rnp - p).normalize_cond();
			PointCloud::Point v = (normal % u).normalize_cond();
			double u_neigh = centered_proj_neigh_point | u;
			double v_neigh = centered_proj_neigh_point | v;
			points2D[i * 2] = u_neigh;
			points2D[i * 2 + 1] = v_neigh;
		}

		GEO::Delaunay_var delaunay = GEO::Delaunay::create(2, "BDEL2d");
		delaunay->set_vertices(ns, points2D);
		PointCloud neighborhood_mesh;
		for (int i = 0; i < proj_neighborhoods[idx].size(); i++)
		{
			PointCloud::Point neigh_point = point_cloud.point(PointCloud::VertexHandle(proj_neighborhoods[idx][i].first));
			neighborhood_mesh.add_vertex(neigh_point);
		}
		for (GEO::index_t t = delaunay->nb_finite_cells(); t < delaunay->nb_cells(); ++t)
		{
			GEO::signed_index_t v0 = delaunay->cell_vertex(t, 0);
			GEO::signed_index_t v1 = delaunay->cell_vertex(t, 1);
			GEO::signed_index_t v2 = delaunay->cell_vertex(t, 2);

			if (v0 != -1 && v1 != -1 && v2 != -1)
			{
				std::vector<PointCloud::VertexHandle>  face_vhandles;
				face_vhandles.push_back(PointCloud::VertexHandle(v0));
				face_vhandles.push_back(PointCloud::VertexHandle(v1));
				face_vhandles.push_back(PointCloud::VertexHandle(v2));
				neighborhood_mesh.add_face(face_vhandles);
			}
		}
		std::vector<PCDNumber> neighborhood_areas;
		compute_mesh_face_areas(neighborhood_mesh, neighborhood_areas);
		PCDNumber vertex_area = compute_mesh_vertex_area(neighborhood_mesh, PointCloud::VertexHandle(0), neighborhood_areas);
		areas[v_it->idx()] = vertex_area;
		PointCloud::Point center_point = neighborhood_mesh.point(PointCloud::VertexHandle(0));
		rings[v_it->idx()].push_back(pair<size_t, PCDNumber>(v_it->idx(), 0));
		PCDNumber max_edge_length = 0;
		for (PointCloud::VertexVertexIter vv_it = neighborhood_mesh.vv_iter(PointCloud::VertexHandle(0)); vv_it.is_valid(); vv_it++)
		{
			PointCloud::Point neigh_point = neighborhood_mesh.point(*vv_it);
			PCDNumber edge_length = (center_point - neigh_point).norm();
			rings[v_it->idx()].push_back(pair<size_t, PCDNumber>(proj_neighborhoods[idx][vv_it->idx()].first, edge_length));
			avg_distance += edge_length;
			edge_count++;
			if (edge_length > max_edge_length)
				max_edge_length = edge_length;
		}
		max_distances_in_rings[v_it->idx()] = max_edge_length;
	}
	avg_distance /= (PCDNumber)(edge_count);


	/*
	****************************************************************************************************************************************************
	Computing regular and big neighborhoods.
	****************************************************************************************************************************************************
	*/
	PCDNumber regular_neighborhood_tolerance = avg_distance * r_r;

	for (size_t i = 0; i < k_filtered_neighborhoods.size(); i++)
	{
		size_t n_elems = 0;
		for (size_t j = 0; j < k_filtered_neighborhoods[i].size(); j++)
		{
			if (k_filtered_neighborhoods[i][j].second <= regular_neighborhood_tolerance)
			{
				regular_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][j].first, k_filtered_neighborhoods[i][j].second));
				n_elems++;
			}
			if (k_filtered_neighborhoods[i][j].second <= (regular_neighborhood_tolerance * r_b))
			{
				big_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][j].first, k_filtered_neighborhoods[i][j].second));
			}
		}
		if (n_elems < min_n_neighbors) //if N_r has a low number of elements, N_r = N_b
		{
			n_elems = 0;
			regular_neighborhoods[i].clear();
			for (size_t j = 0; j < k_filtered_neighborhoods[i].size(); j++)
			{
				if (k_filtered_neighborhoods[i][j].second <= (regular_neighborhood_tolerance*r_b))
				{
					regular_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][j].first, k_filtered_neighborhoods[i][j].second));
					n_elems++;
				}
			}
		}
		if (n_elems < min_n_neighbors) //if N_r and N_b have a low number of elements, N_r = N_20 and N_b = N_20, where N_20 are the first 20 closest points
		{
			// include the first min_n_neighbors
			regular_neighborhoods[i].clear();
			big_neighborhoods[i].clear();
			for (size_t j = 0; j < std::min((int)k_filtered_neighborhoods[i].size(), min_n_neighbors); j++)
			{
				regular_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][j].first, k_filtered_neighborhoods[i][j].second));
				big_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][j].first, k_filtered_neighborhoods[i][j].second));
			}
		}
	}

	/*
	****************************************************************************************************************************************************
	Computing small neighborhoods.
	****************************************************************************************************************************************************
	*/
	for (size_t i = 0; i < k_filtered_neighborhoods.size(); i++)
	{
		size_t n_elems = 0;
		for (size_t j = 0; j < k_filtered_neighborhoods[i].size(); j++)
		{
			if (k_filtered_neighborhoods[i][j].second <= regular_neighborhood_tolerance && k_filtered_neighborhoods[i][j].second <= (max_distances_in_rings[i] * r_s))
			{
				small_neighborhoods[i].push_back(pair<size_t, PCDNumber>(k_filtered_neighborhoods[i][j].first, k_filtered_neighborhoods[i][j].second));
				n_elems++;
			}
		}
		if (n_elems < 5) //small neighborhoods should have at least 5 elements, otherwise, N_s = N_1
		{
			small_neighborhoods[i] = rings[i];
		}
	}

	/*
	****************************************************************************************************************************************************
	Computing regular normals based on regular neighborhoods.
	****************************************************************************************************************************************************
	*/
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		PointCloud::Point t_mean(0, 0, 0);
		PCDMatrix points(regular_neighborhoods[v_it->idx()].size(), 3);
		for (int i = 0; i < regular_neighborhoods[v_it->idx()].size(); i++)
		{
			PointCloud::Point neigh_point = point_cloud.point(PointCloud::VertexHandle(regular_neighborhoods[v_it->idx()][i].first));
			t_mean += neigh_point;
			points(i, 0) = neigh_point[0];
			points(i, 1) = neigh_point[1];
			points(i, 2) = neigh_point[2];
		}
		t_mean /= (double)regular_neighborhoods[v_it->idx()].size();
		for (int i = 0; i < regular_neighborhoods[v_it->idx()].size(); i++)
		{
			points(i, 0) = points(i, 0) - t_mean[0];
			points(i, 1) = points(i, 1) - t_mean[1];
			points(i, 2) = points(i, 2) - t_mean[2];
		}
		PCDMatrix cov = (points.transpose()*points) / (double)k;
		Eigen::JacobiSVD<PCDMatrix> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
		auto U = svd.matrixU();

		PointCloud::Normal normal(U(0, 2), U(1, 2), U(2, 2));
		normal = normal.normalize_cond();
		PointCloud::Normal prevNormal = regular_normals[v_it->idx()];
		prevNormal = prevNormal.normalize_cond();
		PCDNumber dp = normal | prevNormal;
		if (dp < 0.0) //propagation of orientation consistence
		{
			normal = -normal;
		}
		regular_normals[v_it->idx()] = normal;
	}

	/*
	****************************************************************************************************************************************************
	Assignment of regular normals to the input point cloud.
	****************************************************************************************************************************************************
	*/
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		point_cloud.set_normal(*v_it, regular_normals[v_it->idx()]);
	}
}


vector<PointCloud::Normal> normal_filtering(PointCloud & point_cloud, int n_ns, PCDNumber sigma_nn, PCDNumber sigma_ns_ratio, 
	vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals, 
	PCDNumber & avg_distance, bool use_mesh_normals, vector<vector<pair<size_t, PCDNumber>>> & membership_functions)
{
	/*
	****************************************************************************************************************************************************
	In case point_cloud is a mesh, we can use the mesh-based vertex normals, otherwise, use regular normals
	****************************************************************************************************************************************************
	*/
	vector<PointCloud::Normal> normals(regular_normals.size());
	if (use_mesh_normals)
	{
		for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
		{
			PointCloud::Normal normal = (point_cloud.normal(*v_it));
			normals[v_it->idx()] = normal.normalize_cond();
		}
	}
	else
	{
		normals = regular_normals;
	}

	/*
	****************************************************************************************************************************************************
	Computing anisotropic neighborhood normals.
	****************************************************************************************************************************************************
	*/

	PCDNumber tau_u = 0.3;

	vector<PointCloud::Normal> t_normals(membership_functions.size());

	for (int i = 0; i < membership_functions.size(); i++)
	{
		PointCloud::Normal avg_normal = PointCloud::Normal(0, 0, 0);
		vector<int> included_points;
		for (int j = 0; j < membership_functions[i].size(); j++)
		{
			avg_normal += normals[membership_functions[i][j].first] * membership_functions[i][j].second * areas[membership_functions[i][j].first];
			if (membership_functions[i][j].second > tau_u) //thresholding anisotropic neighborhood points
			{
				included_points.push_back(membership_functions[i][j].first);
			}
		}
		avg_normal.normalize_cond();
		t_normals[i] = avg_normal;

		PCDMatrix coord(3, included_points.size());
		for (size_t l = 0; l < included_points.size(); ++l)
		{
			PointCloud::Point inc_point = point_cloud.point(PointCloud::VertexHandle(included_points[l]));
			coord.col(l)[0] = inc_point[0];
			coord.col(l)[1] = inc_point[1];
			coord.col(l)[2] = inc_point[2];
		}

		//computing normal based on thresholded points
		Eigen::Vector3f centroid(coord.row(0).mean(), coord.row(1).mean(), coord.row(2).mean());
		coord.row(0).array() -= centroid(0); coord.row(1).array() -= centroid(1); coord.row(2).array() -= centroid(2);
		auto svd = coord.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::Vector3f plane_normal = svd.matrixU().rightCols<1>();
		PointCloud::Normal new_normal(plane_normal[0], plane_normal[1], plane_normal[2]);
		new_normal = new_normal.normalize_cond();
		if ((new_normal | avg_normal) < 0) //propagates normal orientation consistency
			new_normal = -new_normal;
		t_normals[i] = new_normal;
	}
	normals = t_normals;

	/*
	****************************************************************************************************************************************************
	Bilateral normal smoothing.
	****************************************************************************************************************************************************
	*/

	for (int it = 0; it < n_ns; it++)
	{
		vector<PointCloud::Normal> new_normals = normals;
		for (int i = 0;i < normals.size();i++)
		{
			PCDNumber weight_sum = 0;
			PointCloud::Normal ni = normals[i];
			PointCloud::Point pi = point_cloud.point(PointCloud::VertexHandle(i));
			PointCloud::Normal temp_normal(0, 0, 0);
			for (int j = 0;j < regular_neighborhoods[i].size();j++)
			{
				size_t idx_j = regular_neighborhoods[i][j].first;
				PointCloud::Normal nj = normals[idx_j];
				PointCloud::Point pj = point_cloud.point(PointCloud::VertexHandle(idx_j));
				PCDNumber aj = areas[idx_j];
				PCDNumber spatial_distance = (pi - pj).length();
				PCDNumber spatial_weight = compute_gaussian_weight(spatial_distance, avg_distance*sigma_ns_ratio);
				PCDNumber range_distance = (ni - nj).length();
				PCDNumber range_weight = compute_gaussian_weight(range_distance, sigma_nn);

				PCDNumber weight = spatial_weight * range_weight * aj;
				weight_sum += weight;
				temp_normal += nj * weight;
			}
			temp_normal /= weight_sum;
			temp_normal.normalize_cond();
			new_normals[i] = temp_normal;
		}
		normals = new_normals;
	}

	/*
	****************************************************************************************************************************************************
	Assigning filtered normals to the input point cloud and returing them as vector.
	****************************************************************************************************************************************************
	*/

	vector<PointCloud::Normal> anisotropic_neighborhood_normals = normals;
	for (int i = 0; i < anisotropic_neighborhood_normals.size(); i++)
	{
		point_cloud.set_normal(PointCloud::VertexHandle(i), anisotropic_neighborhood_normals[i]);
	}

	return anisotropic_neighborhood_normals;
}

FEATURE farthest_points(vector<PointCloud::Point> & points, bool use_avg, PCDNumber edge_threshold, vector<PointCloud::Point> & cluster_centroids, 
						vector<int> & cluster_membership)
{
	/*
	****************************************************************************************************************************************************
	Computing mean.
	****************************************************************************************************************************************************
	*/
	PointCloud::Point mean(0, 0, 0);
	for (int i = 0; i < points.size(); i++)
		mean += points[i];
	mean /= (PCDNumber)(points.size());

	/*
	****************************************************************************************************************************************************
	First farthest point computation. Farthest regarding mean.
	****************************************************************************************************************************************************
	*/
	PCDNumber max_dist = MIN_PCDNUMBER;
	vector<int> farthest_points;
	int farthest = 0;
	for (int i = 0; i < points.size(); i++)
	{
		PCDNumber dist = (mean - points[i]).length();
		if (dist > max_dist)
		{
			max_dist = dist;
			farthest = i;
		}
	}
	farthest_points.push_back(farthest);

	/*
	****************************************************************************************************************************************************
	Second farthest point computation. Farthest regarding first farthest point.
	****************************************************************************************************************************************************
	*/
	farthest = 0;
	max_dist = MIN_PCDNUMBER;
	for (int i = 0; i < points.size(); i++)
	{
		PCDNumber dist = (points[farthest_points[0]] - points[i]).length();
		if (dist > max_dist)
		{
			max_dist = dist;
			farthest = i;
		}
	}
	farthest_points.push_back(farthest);
	/*
	****************************************************************************************************************************************************
	Third farthest point computation. Farthest regarding first and second farthest points.
	****************************************************************************************************************************************************
	*/
	farthest = 0;
	max_dist = MIN_PCDNUMBER;
	for (int i = 0; i < points.size(); i++)
	{
		PCDNumber dist = (points[farthest_points[0]] - points[i]).length() + (points[farthest_points[1]] - points[i]).length();
		if (dist > max_dist)
		{
			max_dist = dist;
			farthest = i;
		}
	}
	farthest_points.push_back(farthest);

	/*
	****************************************************************************************************************************************************
	Practical implementation: Fourth farthest point computation. Farthest regarding first, second and third farthest points.
	****************************************************************************************************************************************************
	*/
	farthest = 0;
	max_dist = MIN_PCDNUMBER;
	PCDNumber min_last_farthest = max_dist;
	for (int i = 0; i < points.size(); i++)
	{
		PCDNumber d1 = (points[farthest_points[0]] - points[i]).length();
		PCDNumber d2 = (points[farthest_points[1]] - points[i]).length();
		PCDNumber d3 = (points[farthest_points[2]] - points[i]).length();
		PCDNumber dist = d1 + d2 + d3;
		if (dist > max_dist)
		{
			max_dist = dist;
			farthest = i;
			min_last_farthest = min(min(d1, d2), d3);
		}
	}
	farthest_points.push_back(farthest);

	/*
	****************************************************************************************************************************************************
	Computing maximum edge length and area for the triangle shape in the Gauss sphere. Triangle: (First fp, Second fp, Third fp).
	****************************************************************************************************************************************************
	*/
	PCDNumber max_edge_length = max(max((points[farthest_points[1]] - points[farthest_points[0]]).length(), (points[farthest_points[2]] - points[farthest_points[0]]).length()), (points[farthest_points[0]] - points[farthest_points[2]]).length());

	PointCloud::Point edge1 = points[farthest_points[1]] - points[farthest_points[0]];
	PointCloud::Point edge2 = points[farthest_points[1]] - points[farthest_points[2]];
	PCDNumber area = 0.5f * (edge1 % edge2).length();


	/*
	****************************************************************************************************************************************************
	Computing clusters and feature types
	****************************************************************************************************************************************************
	*/
	FEATURE out_feature_type = FLAT;

	PCDNumber tau_ta = 0.2;

	if (area > tau_ta)
	{
		cluster_membership.clear();
		cluster_centroids.clear();
		cluster_membership.resize(points.size(), 0);
		cluster_centroids.resize(3);
		cluster_centroids[0] = points[farthest_points[0]];
		cluster_centroids[1] = points[farthest_points[1]];
		cluster_centroids[2] = points[farthest_points[2]];
		if (min_last_farthest > 1.0)
		{
			cluster_centroids.push_back(points[farthest_points[3]]);
		}
		for (int i = 0;i < points.size();i++)
		{
			int cluster_id = 0;
			PCDNumber min_d = MAX_PCDNUMBER;
			for (int j = 0;j < cluster_centroids.size();j++)
			{
				PCDNumber d = (cluster_centroids[j] - points[i]).norm();
				if (d < min_d)
				{
					cluster_id = j;
					min_d = d;
				}
			}
			cluster_membership[i] = cluster_id;
		}
		out_feature_type = CORNER; //corner
	}
	else if (max_edge_length > edge_threshold)
	{
		cluster_membership.clear();
		cluster_centroids.clear();
		cluster_membership.resize(points.size(), 0);
		cluster_centroids.resize(2);
		cluster_centroids[0] = points[farthest_points[0]];
		cluster_centroids[1] = points[farthest_points[1]];
		for (int i = 0;i < points.size();i++)
		{
			int cluster_id = 0;
			PCDNumber min_d = MAX_PCDNUMBER;
			for (int j = 0;j < cluster_centroids.size();j++)
			{
				PCDNumber d = (cluster_centroids[j] - points[i]).norm();
				if (d < min_d)
				{
					cluster_id = j;
					min_d = d;
				}
			}
			cluster_membership[i] = cluster_id;
		}
		out_feature_type = EDGE; //edge
	}
	else
	{
		cluster_centroids.clear();
		cluster_membership.clear();
		cluster_centroids.push_back(points[farthest_points[0]]);
		cluster_membership.resize(points.size(), 0);
		out_feature_type = FLAT; //flat
	}

	/*
	****************************************************************************************************************************************************
	Use centroids instead of seeds (initial farthest points)
	****************************************************************************************************************************************************
	*/
	if (use_avg)
	{
		vector<PointCloud::Point> avg_centroids(cluster_centroids.size(), PointCloud::Point(0, 0, 0));
		vector<double> avg_centroids_count(cluster_centroids.size(), 0.0);
		for (size_t i = 0;i < points.size();i++)
		{
			size_t cluster_idx = cluster_membership[i];
			avg_centroids[cluster_idx] += points[i];
			avg_centroids_count[cluster_idx] += 1.0;
		}
		for (size_t i = 0;i < cluster_centroids.size();i++)
		{
			avg_centroids[i] /= avg_centroids_count[i];
			cluster_centroids[i] = avg_centroids[i];
		}
	}

	return out_feature_type;
}

void normal_correction(PointCloud & point_cloud, vector<PointCloud::Normal> &normals, vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods, 
	vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals, 
	PCDNumber & avg_distance, GEO::NearestNeighborSearch_var & NN, PCDNumber tau_n, int n_nc, vector<vector<PointCloud::Normal>> & multi_normals)
{
	vector<FEATURE> classification(normals.size(), FLAT);

	multi_normals.clear();
	multi_normals.resize(normals.size());

	for (int i = 0;i < normals.size();i++)
	{
		multi_normals[i].push_back(normals[i]);
	}

	size_t n_normal_correction_passes = n_nc;

	PCDNumber e_mn_ratio = 0.1;
	/*
	****************************************************************************************************************************************************
	Normal correction iterations
	****************************************************************************************************************************************************
	*/

	for (size_t pass = 0;pass < n_normal_correction_passes;pass++)
	{
		vector<PointCloud::Normal> corrected_normals = normals;

		for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
		{
			PointCloud::Point p = point_cloud.point(*v_it);
			PointCloud::Normal normal = normals[v_it->idx()].normalize_cond();
			/*
			****************************************************************************************************************************************************
			Selecting data for clustering (without including the evaluated point)
			****************************************************************************************************************************************************
			*/
			vector<PointCloud::Point> neigh_normals;
			vector<PointCloud::Point> neigh_points;
			vector<size_t> neigh_ids;
			for (int i = 1; i < regular_neighborhoods[v_it->idx()].size(); i++) // ignore the current vertex
			{
				int id_i = regular_neighborhoods[v_it->idx()][i].first;
				PointCloud::Normal n_i = normals[id_i].normalize_cond();
				PointCloud::Point p_i = point_cloud.point(PointCloud::VertexHandle(id_i));
				neigh_normals.push_back(n_i);
				neigh_points.push_back(p_i);
				neigh_ids.push_back(id_i);
			}
			/*
			****************************************************************************************************************************************************
			Rough feature classification (candidate selection)
			****************************************************************************************************************************************************
			*/
			PCDNumber max_normal_difference = 0;
			for (int i = 0; i < small_neighborhoods[v_it->idx()].size(); i++)
			{
				PointCloud::Normal n_i = normals[small_neighborhoods[v_it->idx()][i].first].normalize_cond();
				for (int j = i + 1; j < small_neighborhoods[v_it->idx()].size(); j++)
				{
					PointCloud::Normal n_j = normals[small_neighborhoods[v_it->idx()][j].first].normalize_cond();
					PCDNumber normal_difference = (n_j - n_i).norm();
					if (normal_difference > max_normal_difference)
						max_normal_difference = normal_difference;
				}
			}

			if ((max_normal_difference) > tau_n)
			{
				classification[v_it->idx()] = EDGE;
			}
			else
				classification[v_it->idx()] = FLAT;

			/*
			****************************************************************************************************************************************************
			Clustering based on Gauss sphere farthest points estimation, normal correction, and multi-normal computation
			****************************************************************************************************************************************************
			*/
			vector<PointCloud::Point> cluster_centroids;
			vector<int> cluster_membership;
			//Defining clusters
			if (classification[v_it->idx()] == EDGE)
			{
				classification[v_it->idx()] = farthest_points(neigh_normals, true, tau_n, cluster_centroids, cluster_membership);
			}
			//Computing cluster elements, cluster avg location, avg normal and avg distance regarding per-point planes
			if (classification[v_it->idx()] == EDGE || classification[v_it->idx()] == CORNER)
			{
				vector<PCDNumber> cluster_avg_distances(cluster_centroids.size(), 0);
				vector<PointCloud::Point> cluster_location_centroids(cluster_centroids.size(), PointCloud::Point(0, 0, 0));
				vector<PointCloud::Normal> cluster_normal_centroids(cluster_centroids.size(), PointCloud::Normal(0, 0, 0));
				vector<PCDNumber> cluster_count(cluster_centroids.size(), 0);
				vector<vector<PointCloud::Point>> cluster_points(cluster_centroids.size());
				for (int i = 0;i < neigh_points.size();i++)
				{
					PCDNumber d = point_to_plane_distance(p, neigh_normals[i], neigh_points[i]);
					cluster_avg_distances[cluster_membership[i]] += d;
					cluster_location_centroids[cluster_membership[i]] += neigh_points[i];
					cluster_normal_centroids[cluster_membership[i]] += neigh_normals[i];
					cluster_count[cluster_membership[i]] += 1.0;
					cluster_points[cluster_membership[i]].push_back(neigh_points[i]);
				}
				for (int i = 0;i < cluster_avg_distances.size();i++)
				{
					cluster_avg_distances[i] = cluster_avg_distances[i] / cluster_count[i];
					cluster_location_centroids[i] /= cluster_count[i];
					cluster_normal_centroids[i] /= cluster_count[i];
				}
				PCDNumber min_d_cluster = MAX_PCDNUMBER;
				PCDNumber min_d_plane = MAX_PCDNUMBER;
				int min_d_cluster_idx = 0;
				vector<PCDNumber> t_distances;
				//Computing minimum d_cluster and assigning multi-normals based on d_plane only
				for (int i = 0;i < cluster_avg_distances.size();i++)
				{
					
					PCDNumber d_plane = cluster_avg_distances[i];
					PCDNumber d_ch = point_to_plane_region_distance(p, cluster_normal_centroids[i], cluster_location_centroids[i], cluster_points[i]);
					PCDNumber d_cluster = d_plane + d_ch;
					t_distances.push_back(d_plane);
					if (d_plane < min_d_plane)
					{
						min_d_plane = d_plane;
					}

					if (d_cluster < min_d_cluster)
					{
						min_d_cluster = d_cluster;
						min_d_cluster_idx = i;
					}
				}
				//Correcting normals 
				PointCloud::Normal min_d_cluster_normal = cluster_centroids[min_d_cluster_idx];
				if ((min_d_cluster_normal - normal).norm() > tau_n)
				{
					corrected_normals[v_it->idx()] = min_d_cluster_normal;
				}

				//Assigning multi-normals
				multi_normals[v_it->idx()].clear();
				vector<size_t> indexes = sort_indexes(t_distances);
				PCDNumber distance_tolerance = min_d_plane + avg_distance * e_mn_ratio;
				for (int i = 0;i < indexes.size();i++)
				{
					if (t_distances[indexes[i]] < distance_tolerance)
						multi_normals[v_it->idx()].push_back(cluster_centroids[indexes[i]]);
				}
			}
		}
		/*
		****************************************************************************************************************************************************
		Updating normals vector and input point cloud normals using corrected normals
		****************************************************************************************************************************************************
		*/
		for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
		{
			normals[v_it->idx()] = corrected_normals[v_it->idx()];
			point_cloud.set_normal(*v_it, corrected_normals[v_it->idx()]);
		}
	}
}




void classify_points(PCDNumber tau_n, PCDNumber theta_half, PointCloud & point_cloud, vector<PointCloud::Normal> & normals, 
	vector<vector<pair<size_t, PCDNumber>>> & regular_neighborhoods, vector<vector<pair<size_t, PCDNumber>>> & big_neighborhoods, 
	vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals,
	PCDNumber & avg_distance, GEO::NearestNeighborSearch_var & NN, PCDNumber delta_cc_ratio, vector<FEATURE> & classification, 
	vector<PointCloud::Normal> & feature_dirs)
{
	classification.clear();
	classification.resize(normals.size(), FLAT);
	feature_dirs.clear();
	feature_dirs.resize(normals.size());

	vector<int> candidate_feature_points;

	candidate_feature_points.clear();
	std::fill(classification.begin(), classification.end(), FLAT);

	vector<vector<int>> all_cluster_memberships(point_cloud.n_vertices());
	vector<vector<vector<int>>> all_cluster_memberships_ids(point_cloud.n_vertices());
	vector<vector<PointCloud::Point>> all_cluster_centroids(point_cloud.n_vertices());
	vector<vector<PointCloud::Point>> all_cluster_location_centroids(point_cloud.n_vertices());

	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		PointCloud::Normal normal = normals[v_it->idx()].normalize_cond();
		/*
		****************************************************************************************************************************************************
		Selecting data for clustering (including the evaluated point)
		****************************************************************************************************************************************************
		*/
		vector<PointCloud::Point> neigh_normals;
		vector<PointCloud::Point> neigh_points;
		for (int i = 0; i < small_neighborhoods[v_it->idx()].size(); i++)
		{
			PointCloud::Normal n_i = normals[small_neighborhoods[v_it->idx()][i].first].normalize_cond();
			PointCloud::Point p_i = point_cloud.point(PointCloud::VertexHandle(small_neighborhoods[v_it->idx()][i].first));
			neigh_normals.push_back(n_i);
			neigh_points.push_back(p_i);
		}

		/*
		****************************************************************************************************************************************************
		Rough feature classification (candidate selection)
		****************************************************************************************************************************************************
		*/
		PCDNumber max_normal_difference = 0;

		for (int i = 0; i < small_neighborhoods[v_it->idx()].size(); i++)
		{
			PointCloud::Normal n_i = normals[small_neighborhoods[v_it->idx()][i].first].normalize_cond();
			for (int j = i + 1; j < small_neighborhoods[v_it->idx()].size(); j++)
			{
				PointCloud::Normal n_j = normals[small_neighborhoods[v_it->idx()][j].first].normalize_cond();
				PCDNumber normal_difference = (n_j - n_i).norm();
				if (normal_difference > max_normal_difference)
					max_normal_difference = normal_difference;
			}
		}

		if ((max_normal_difference) > tau_n)
		{
			classification[v_it->idx()] = EDGE;
		}
		else
			classification[v_it->idx()] = FLAT;

		/*
		****************************************************************************************************************************************************
		Clustering based on Gauss sphere farthest points estimation, normal correction, and multi-normal computation
		****************************************************************************************************************************************************
		*/
		vector<PointCloud::Point> cluster_centroids;
		vector<int> cluster_membership;
		if (classification[v_it->idx()] == EDGE)
		{
			classification[v_it->idx()] = farthest_points(neigh_normals, true, tau_n, cluster_centroids, cluster_membership);

			if (classification[v_it->idx()] != FLAT)
				candidate_feature_points.push_back(v_it->idx());
			all_cluster_memberships[v_it->idx()] = cluster_membership;
			all_cluster_centroids[v_it->idx()] = cluster_centroids;

			vector<vector<int>> per_cluster_ids(cluster_centroids.size());
			vector<vector<int>> cluster_spatial_centroids(cluster_centroids.size());
			for (size_t k = 0;k < cluster_membership.size();k++)
			{
				per_cluster_ids[cluster_membership[k]].push_back(small_neighborhoods[v_it->idx()][k].first);
			}
			all_cluster_memberships_ids[v_it->idx()] = per_cluster_ids;

			vector<PointCloud::Point> cluster_localization_centroids(cluster_centroids.size());
			for (size_t k = 0; k < per_cluster_ids.size();k++)
			{
				PointCloud::Point localization_centroid(0, 0, 0);
				PCDNumber loc_cen_count = 0;
				for (size_t l = 0; l < per_cluster_ids[k].size(); l++)
				{
					localization_centroid += point_cloud.point(PointCloud::VertexHandle(per_cluster_ids[k][l]));
					loc_cen_count += 1.0;
				}
				localization_centroid /= loc_cen_count;
				cluster_localization_centroids[k] = localization_centroid;
			}
			all_cluster_location_centroids[v_it->idx()] = cluster_localization_centroids;
		}
	}

	/*
	****************************************************************************************************************************************************
	Creating a smooth version of the point cloud for convexity analysis (Laplacian smoothing)
	****************************************************************************************************************************************************
	*/
	int n_smoothing_iterations = 10;
	PCDNumber smoothing_step_size = 0.2;

	PointCloud smooth_point_cloud_cc = point_cloud;
	for (int it = 0;it < n_smoothing_iterations;it++)
	{
		vector<PointCloud::Point> new_positions(smooth_point_cloud_cc.n_vertices());
		for (PointCloud::VertexIter v_it = smooth_point_cloud_cc.vertices_begin(); v_it != smooth_point_cloud_cc.vertices_end(); v_it++)
		{
			PointCloud::Point p = smooth_point_cloud_cc.point(*v_it);
			size_t idx = v_it->idx();
			PointCloud::Point centroid(0, 0, 0);
			PCDNumber count = 0;
			for (int i = 1;i < small_neighborhoods[idx].size();i++)
			{
				size_t neigh_idx = small_neighborhoods[idx][i].first;
				centroid += smooth_point_cloud_cc.point(PointCloud::VertexHandle(neigh_idx));
				count += 1.0;
			}
			centroid /= count;
			new_positions[idx] = p + smoothing_step_size *(centroid - p);
		}
		for (size_t i = 0;i < new_positions.size();i++)
		{
			smooth_point_cloud_cc.set_point(PointCloud::VertexHandle(i), new_positions[i]);
		}
	}


	/*
	****************************************************************************************************************************************************
	Convexity analysis
	****************************************************************************************************************************************************
	*/
	PCDNumber e_cc_ratio = 0.2;
	vector<int> concave_undefined_convex(smooth_point_cloud_cc.n_vertices(), 0); //concave: -1, undefined: 0, convex: 1
	for (PointCloud::VertexIter v_it = smooth_point_cloud_cc.vertices_begin(); v_it != smooth_point_cloud_cc.vertices_end(); v_it++)
	{
		/*
		****************************************************************************************************************************************************
		Computing a smooth normal
		****************************************************************************************************************************************************
		*/
		PointCloud::Point p = smooth_point_cloud_cc.point(*v_it);
		PointCloud::Normal regular_normal = regular_normals[v_it->idx()];
		PointCloud::Normal avg_regular_normal(0, 0, 0);
		for (int i = 0; i < small_neighborhoods[v_it->idx()].size(); i++)
		{
			size_t neigh_idx = small_neighborhoods[v_it->idx()][i].first;
			avg_regular_normal += regular_normals[neigh_idx];
		}
		avg_regular_normal = avg_regular_normal.normalize_cond();
		/*
		****************************************************************************************************************************************************
		If point is FLAT, set as undefined. Otherwise, continue analysis.
		****************************************************************************************************************************************************
		*/
		if (classification[v_it->idx()] == CORNER || classification[v_it->idx()] == EDGE)
		{
			/*
			****************************************************************************************************************************************************
			Computing big neighborhood centroid (uniform)
			****************************************************************************************************************************************************
			*/
			PointCloud::Point neighborhood_centroid(0, 0, 0);
			double sum_dist = 0;
			double sum_weight = 0;
			for (int i = 1; i < big_neighborhoods[v_it->idx()].size(); i++)
			{
				size_t neigh_idx = (big_neighborhoods[v_it->idx()][i]).first;
				PointCloud::Point neigh_point = smooth_point_cloud_cc.point(PointCloud::VertexHandle(neigh_idx));
				PCDNumber neigh_dist = (neigh_point - p).length();
				neighborhood_centroid += neigh_point;
				sum_weight += 1.0;
			}
			neighborhood_centroid = neighborhood_centroid / sum_weight;

			/*
			****************************************************************************************************************************************************
			Measuring convexity/concavity using the signed distance from the evaluated point to the plane described by the smooth normal and the centroid
			****************************************************************************************************************************************************
			*/
			PCDNumber l = point_to_plane_distance(p, avg_regular_normal, neighborhood_centroid);
			PCDNumber e_cc = avg_distance * e_cc_ratio;
			PCDNumber sign = (avg_regular_normal | (p - neighborhood_centroid).normalize_cond());
			if (sign > 0 && l > e_cc) //convex
			{
				concave_undefined_convex[v_it->idx()] = 1;
			}
			else if (sign < 0 && l>e_cc) //concave
			{
				concave_undefined_convex[v_it->idx()] = -1;
			}
			else //undefined
			{
				concave_undefined_convex[v_it->idx()] = 0;
			}
		}
		else //undefined
		{
			concave_undefined_convex[v_it->idx()] = 0;
		}
	}


	/*
	****************************************************************************************************************************************************
	Feature classification and estimation of feature lines (edge directions)
	****************************************************************************************************************************************************
	*/
	std::fill(classification.begin(), classification.end(), FLAT);
	PCDNumber delta_cc = avg_distance * delta_cc_ratio;
	for (size_t i = 0; i < candidate_feature_points.size(); i++) //evaluating just candidate feature points
	{
		/*
		****************************************************************************************************************************************************
		Computing displacement amounts for edge lines estimation, based on convexity analysis
		****************************************************************************************************************************************************
		*/
		size_t idx1 = candidate_feature_points[i];
		PCDNumber delta_cc_local = delta_cc;
		if (concave_undefined_convex[idx1] == -1)
		{
			delta_cc_local *= -1.0;
		}
		else if (concave_undefined_convex[idx1] == 0)
		{
			delta_cc_local *= 0.0;
		}

		/*
		****************************************************************************************************************************************************
		Initializing data for feature classification
		****************************************************************************************************************************************************
		*/
		PointCloud::Normal n1 = normals[idx1];
		PointCloud::Normal reg_n1 = regular_normals[idx1];
		PointCloud::Point p1 = point_cloud.point(PointCloud::VertexHandle(idx1));

		int cluster_id_i = all_cluster_memberships[idx1][0];
		vector<int> same_cloud_ids = all_cluster_memberships_ids[idx1][cluster_id_i];
		PointCloud::Point p1_cluster = all_cluster_location_centroids[idx1][cluster_id_i];
		int global_feature_count = 0;

		PointCloud::Point opt_feature_dir(0, 0, 0);
		double min_dist_opt_feature_dir = MAX_PCDNUMBER;
		/*
		****************************************************************************************************************************************************
		Per-cluster analysis
		****************************************************************************************************************************************************
		*/
		for (size_t k = 0;k < all_cluster_memberships_ids[idx1].size();k++)
		{
			double feature_count = 0.0;
			double total_count = 0.0;
			PointCloud::Point feature_dir(0, 0, 0);

			if (k == cluster_id_i) continue; //just evaluate the other clusters

			/*
			****************************************************************************************************************************************************
			Per-point in cluster analysis
			****************************************************************************************************************************************************
			*/
			for (size_t j = 0; j < all_cluster_memberships_ids[idx1][k].size();j++)
			{
				size_t idx2 = all_cluster_memberships_ids[idx1][k][j];
				PointCloud::Normal n2 = normals[idx2];
				PointCloud::Normal reg_n2 = regular_normals[idx2];
				PointCloud::Point p2 = point_cloud.point(PointCloud::VertexHandle(idx2));
				PointCloud::Point p2_cluster = all_cluster_location_centroids[idx1][k];
				PointCloud::Point lv;
				PointCloud::Point lp;
				plane_to_plane_intersection(p1 + n1 * delta_cc_local, n1, p2 + n2 * delta_cc_local, n2, lp, lv);
				lv = lv.normalize_cond();
				PointCloud::Point closest_point_1;
				PCDNumber intersection_distance_1 = point_to_line(p1, lp, lv, closest_point_1); //distance from p_i to the edge line e_ij
				PCDNumber intersection_distance_2 = (p2 - closest_point_1).norm(); //distance from p_j to the edge line e_ij
				PCDNumber agd1 = angular_distance_on_tangent_plane(n1, lv, p1, n1, p2); //angular distance on plane P_i
				PCDNumber agd2 = angular_distance_on_tangent_plane(n2, lv, p2, n2, p2); //angular distance on plane P_j
				if (!(agd1 < theta_half || agd2 < theta_half)) 
				{
					/*If the other cluster point is not in the narrow neighborhood, just ignore it*/
					continue;
				}

				/*Computing the closest point to e_ij within the evaluated point cluster C_1*/
				PCDNumber min_d_l = intersection_distance_1;
				PCDNumber min_d_l_idx = idx1;
				for (size_t l = 0;l < same_cloud_ids.size();l++)
				{
					PointCloud::Point closest_point_l;
					size_t l_idx = same_cloud_ids[l];
					if (l_idx == idx1) continue;
					PointCloud::Point l_point = point_cloud.point(PointCloud::VertexHandle(l_idx));

					if (angular_distance_on_tangent_plane(n1, lv, p1, n1, l_point) > theta_half)
					{
						/*If the same cluster point is not in the narrow neighborhood, just ignore it*/
						continue;
					}
					PCDNumber d_l = (closest_point_1 - l_point).norm();
					if (d_l < min_d_l)
					{
						min_d_l = d_l;
						min_d_l_idx = l_idx;
					}
				}

				/*Check if p_i is the closest point to e_ij within the narrow neighborhood*/
				if (intersection_distance_1 < intersection_distance_2 && min_d_l_idx == idx1)
				{
					feature_count += 1.0;
					feature_dir += PointCloud::Point(lv[0], lv[1], lv[2]);
					if (intersection_distance_1 < min_dist_opt_feature_dir)
					{
						/*Considering the closest e_ij as opt_feature_dir*/
						min_dist_opt_feature_dir = intersection_distance_1;
						opt_feature_dir = PointCloud::Point(lv[0], lv[1], lv[2]);
					}
				}
				total_count += 1.0;
			}
			/*
			****************************************************************************************************************************************************
			Classifying the evaluated point as edge if it is the closest point regarding all the points included in the narrow neighborhood 
			defined by a single cluster. Feature count denotes the number of times the evaluated point p_i was considered as the closest one
			regarding a point p_j in another cluster. Total count denotes the number of points within the other cluster used for analysis.
			****************************************************************************************************************************************************
			*/
			double ratio = 0.0;
			if (total_count != 0)
				ratio = feature_count / total_count;

			feature_dir = feature_dir.normalize_cond();
			if (ratio > 0.98)
			{
				classification[idx1] = EDGE;
				global_feature_count++;
			}
		}
		opt_feature_dir = opt_feature_dir.normalize_cond();

		if (global_feature_count >= 2)
		{
			classification[idx1] = CORNER;
			bool is_maximum = true;
			bool is_minimum = true;
			for (int j = 1; j < big_neighborhoods[idx1].size(); j++)
			{
				size_t neigh_idx = (big_neighborhoods[idx1][j]).first;
				PointCloud::Point neigh_point = point_cloud.point(PointCloud::VertexHandle(neigh_idx));
				double dist = reg_n1 | (neigh_point - p1);
				if (dist > 0)
				{
					is_maximum = false;
				}
				if (dist < 0)
				{
					is_minimum = false;
				}
			}
			if (is_maximum || is_minimum)
			{
				classification[idx1] = CORNER;
			}
			else
			{
				classification[idx1] = EDGE;
			}
		}
		feature_dirs[idx1] = opt_feature_dir;
	}
}

void update_point_positions(PointCloud & point_cloud, PointCloud & original_point_cloud, vector<PointCloud::Normal> & normals, 
	vector<vector<pair<size_t, PCDNumber>>> & small_neighborhoods, vector<PCDNumber> & areas, vector<PointCloud::Normal> & regular_normals, 
	PCDNumber & avg_distance, GEO::NearestNeighborSearch_var & NN, int n_fp, PCDNumber upsilon_f, 
	PCDNumber upsilon_ec, PCDNumber tau_o_ratio, vector<FEATURE> & classification, vector<PointCloud::Normal> & feature_dirs,
	vector<vector<PointCloud::Normal>> & multi_normals)
{
	PCDNumber sigma_ps_ratio = 2.0;
	PCDNumber sigma_pn = 0.5;

	PCDNumber tolerance = avg_distance * tau_o_ratio;

	vector<PointCloud::Point> new_positions(point_cloud.n_vertices());

	/*
	****************************************************************************************************************************************************
	Flat points update
	****************************************************************************************************************************************************
	*/
	for (int it = 0;it < n_fp;it++)
	{
		for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
		{
			PointCloud::Point p = point_cloud.point(*v_it);
			PointCloud::Point op = original_point_cloud.point(*v_it);
			PointCloud::Normal normal = normals[v_it->idx()];

			if (classification[v_it->idx()] == FLAT) // FLAT
			{
				PointCloud::Point acc(0, 0, 0);
				double accweight = 0;
				PointCloud::Normal np = normals[v_it->idx()].normalize_cond();
				for (int i = 1; i < small_neighborhoods[v_it->idx()].size(); i++)
				{
					int idq = small_neighborhoods[v_it->idx()][i].first;
					PointCloud::Point q = point_cloud.point(PointCloud::VertexHandle(small_neighborhoods[v_it->idx()][i].first));
					PointCloud::Normal nq = normals[small_neighborhoods[v_it->idx()][i].first].normalize_cond();
					double spatial_weight = compute_gaussian_weight((q - p).norm(), avg_distance*sigma_ps_ratio);
					double normal_weight = compute_gaussian_weight((nq - np).norm(), 0.5);
					for (int mni = 0;mni < multi_normals[idq].size();mni++) // remove this for to ignore multi normal
					{
						PointCloud::Normal t_nq = multi_normals[idq][mni];
						double t_normal_weight = compute_gaussian_weight((t_nq - np).norm(), 0.5);
						if (t_normal_weight > normal_weight)
						{
							normal_weight = t_normal_weight;
							nq = t_nq;
						}
					}
					double weight = normal_weight * spatial_weight;
					double dot_prod = nq | (q - p);
					acc += weight * (dot_prod)*np;
					accweight += weight;
				}
				acc = acc / accweight;
				PointCloud::Point virtual_position = p + 1.0*acc;
				PointCloud::Point displacement = virtual_position - p;
				PointCloud::Point new_position = p + displacement * upsilon_f;
				if ((new_position - op).norm() < tolerance)
				{
					new_positions[v_it->idx()] = new_position;
				}
				else
				{
					new_positions[v_it->idx()] = p;
				}
			}
			else
			{
				new_positions[v_it->idx()] = p;
			}
		}
		for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
		{
			point_cloud.set_point(*v_it, new_positions[v_it->idx()]);
		}
	}

	/*
	****************************************************************************************************************************************************
	Corner and edge points update
	****************************************************************************************************************************************************
	*/
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		PointCloud::Point op = original_point_cloud.point(*v_it);
		PointCloud::Normal normal = normals[v_it->idx()];
		if (classification[v_it->idx()] == CORNER)
		{
			PCDMatrix t1 = PCDMatrix::Zero(3, 3);
			PCDMatrix t2 = PCDMatrix::Zero(3, 3);
			PCDMatrix normal_matrix(3, 1);
			normal_matrix << normal[0], normal[1], normal[2];
			for (int i = 1; i < small_neighborhoods[v_it->idx()].size(); i++)
			{
				size_t neigh_idx = (small_neighborhoods[v_it->idx()][i]).first;
				PCDMatrix temp = PCDMatrix::Zero(3, 3);
				PointCloud::Normal neigh_normal = normals[neigh_idx];
				PCDMatrix neigh_normal_matrix(3, 1);
				neigh_normal_matrix << neigh_normal[0], neigh_normal[1], neigh_normal[2];
				PointCloud::Point neigh_point = point_cloud.point(PointCloud::VertexHandle(neigh_idx));
				PCDMatrix neigh_point_matrix(3, 1);
				neigh_point_matrix << neigh_point[0], neigh_point[1], neigh_point[2];
				temp = neigh_normal_matrix * neigh_normal_matrix.transpose();
				t1 += temp;
				temp *= neigh_point_matrix;
				t2 += temp;
			}
			PCDMatrix new_position_matrix = t1.inverse() * t2;
			PointCloud::Point virtual_position(new_position_matrix(0, 0), new_position_matrix(1, 0), new_position_matrix(2, 0));
			PointCloud::Point displacement = virtual_position - p;
			PointCloud::Point new_position = p + displacement * upsilon_ec;
			if ((new_position - op).norm() < tolerance)
			{
				new_positions[v_it->idx()] = new_position;
			}
			else
			{
				new_positions[v_it->idx()] = p;
			}
		}
		else if (classification[v_it->idx()] == EDGE)
		{
			PointCloud::Normal feature_dir = feature_dirs[v_it->idx()];
			PCDMatrix feature_dir_matrix(3, 1);
			feature_dir_matrix << feature_dir[0], feature_dir[1], feature_dir[2];
			PCDMatrix point_matrix(3, 1);
			point_matrix << p[0], p[1], p[2];
			PCDMatrix t1 = PCDMatrix::Zero(3, 3);
			PCDMatrix t2 = PCDMatrix::Zero(3, 3);
			for (int i = 1; i < small_neighborhoods[v_it->idx()].size(); i++)
			{
				size_t neigh_idx = (small_neighborhoods[v_it->idx()][i]).first;
				PointCloud::Point neigh_point = point_cloud.point(PointCloud::VertexHandle(neigh_idx));
				PCDMatrix neigh_point_matrix(3, 1);
				neigh_point_matrix << neigh_point[0], neigh_point[1], neigh_point[2];
				PointCloud::Normal neigh_normal = normals[neigh_idx];
				PointCloud::Normal proj_neigh_normal = neigh_normal - (neigh_normal | feature_dir) * feature_dir;
				PCDMatrix proj_neigh_normal_matrix(3, 1);
				proj_neigh_normal_matrix << proj_neigh_normal[0], proj_neigh_normal[1], proj_neigh_normal[2];
				t1 += proj_neigh_normal_matrix * proj_neigh_normal_matrix.transpose() + feature_dir_matrix * feature_dir_matrix.transpose();
				t2 += proj_neigh_normal_matrix * proj_neigh_normal_matrix.transpose() * neigh_point_matrix + feature_dir_matrix * feature_dir_matrix.transpose() * point_matrix;
			}
			PCDMatrix new_position_matrix = t1.inverse() * t2;
			PointCloud::Point virtual_position(new_position_matrix(0, 0), new_position_matrix(1, 0), new_position_matrix(2, 0));
			PointCloud::Point displacement = virtual_position - p;
			PointCloud::Point new_position = p + displacement * upsilon_ec;
			if ((new_position - op).norm() < tolerance)
			{
				new_positions[v_it->idx()] = new_position;
			}
			else
			{
				new_positions[v_it->idx()] = p;
			}
		}
		else // FLAT
		{
			new_positions[v_it->idx()] = p;
		}

	}
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		point_cloud.set_point(*v_it, new_positions[v_it->idx()]);
	}
}

PointCloud denoise(PointCloud & point_cloud, DenoisingParameters & denoising_parameters, DenoisingData & denoising_data)
{

	PointCloud processed_point_cloud = point_cloud;
	PCDNumber scale_factor;
	normalize_point_cloud(processed_point_cloud, scale_factor);
	rescale_point_cloud(point_cloud,scale_factor);

	PCDNumber prepare_time = 0.0;
	PCDNumber anisotropic_neighborhoods_time = 0.0;
	PCDNumber normal_filtering_time = 0.0;
	PCDNumber normal_correction_time = 0.0;
	PCDNumber feature_classification_time = 0.0;
	PCDNumber point_update_time = 0.0;



	clock_t begin = clock();
	clock_t t_begin;
	clock_t t_end;

	for (int ext_it = 0;ext_it<denoising_parameters.n_ext;ext_it++)
	{
		cout << "External iteration number: " << ext_it << endl;
		cout << "\tPrepare point cloud" << endl;
		t_begin = clock();
		prepare_point_cloud_data(processed_point_cloud, denoising_parameters.k, denoising_parameters.r_r, denoising_data.regular_neighborhoods,
			denoising_data.big_neighborhoods, denoising_data.small_neighborhoods, denoising_data.areas, denoising_data.regular_normals,
			denoising_data.avg_distance, denoising_data.NN);
		t_end = clock();
		prepare_time += PCDNumber(t_end - t_begin) / CLOCKS_PER_SEC;

		cout << "\tAnisotropic neighborhoods" << endl;
		t_begin = clock();
		compute_all_anisotropic_neighborhoods(processed_point_cloud, denoising_data.areas, denoising_data.regular_normals,
			denoising_data.regular_neighborhoods, denoising_data.small_neighborhoods, denoising_parameters.alpha, denoising_parameters.beta,
			denoising_parameters.gamma, denoising_parameters.a_0_ratio, denoising_data.adaptive_kernels);
		t_end = clock();
		anisotropic_neighborhoods_time += PCDNumber(t_end - t_begin) / CLOCKS_PER_SEC;

		for (int int_it = 0; int_it < denoising_parameters.n_int; int_it++)
		{
			cout << "\tInternal iteration number: " << int_it << endl;
			cout << "\t\tNormal filtering" << endl;
			t_begin = clock();
			denoising_data.normals = normal_filtering(processed_point_cloud, denoising_parameters.n_ns, denoising_parameters.sigma_nn,
				denoising_parameters.sigma_ns_ratio, denoising_data.regular_neighborhoods, denoising_data.areas, denoising_data.regular_normals, 
				denoising_data.avg_distance, false, denoising_data.adaptive_kernels);
			t_end = clock();
			normal_filtering_time += PCDNumber(t_end - t_begin) / CLOCKS_PER_SEC;

			cout << "\t\tNormal correction" << endl;
			t_begin = clock();
			normal_correction(processed_point_cloud, denoising_data.normals, denoising_data.regular_neighborhoods,
				denoising_data.small_neighborhoods, denoising_data.areas, denoising_data.regular_normals, 
				denoising_data.avg_distance, denoising_data.NN, denoising_parameters.tau_n, denoising_parameters.n_nc, denoising_data.multi_normals);
			t_end = clock();
			normal_correction_time += PCDNumber(t_end - t_begin) / CLOCKS_PER_SEC;

			cout << "\t\tFeature classification" << endl;
			t_begin = clock();
			classify_points(denoising_parameters.tau_n, denoising_parameters.theta_half, processed_point_cloud, denoising_data.normals, 
				denoising_data.regular_neighborhoods, denoising_data.big_neighborhoods, denoising_data.small_neighborhoods, 
				denoising_data.areas, denoising_data.regular_normals, denoising_data.avg_distance, denoising_data.NN,
				denoising_parameters.delta_cc_ratio, denoising_data.feature_classification, denoising_data.feature_dirs);
			t_end = clock();
			feature_classification_time += PCDNumber(t_end - t_begin) / CLOCKS_PER_SEC;

			cout << "\t\tPoint updating" << endl;
			t_begin = clock();
			update_point_positions(processed_point_cloud, point_cloud,
				denoising_data.normals,
				denoising_data.small_neighborhoods, denoising_data.areas, denoising_data.regular_normals,
				denoising_data.avg_distance, denoising_data.NN, denoising_parameters.n_fp,
				denoising_parameters.upsilon_f, denoising_parameters.upsilon_ec,
				denoising_parameters.tau_o_ratio, denoising_data.feature_classification, denoising_data.feature_dirs,
				denoising_data.multi_normals);
			t_end = clock();
			point_update_time += PCDNumber(t_end - t_begin) / CLOCKS_PER_SEC;
		}
	}

	retrieve_point_cloud(processed_point_cloud, scale_factor);
	retrieve_point_cloud(point_cloud, scale_factor);

	int int_n_ops = denoising_parameters.n_int * denoising_parameters.n_ext;
	int ext_n_ops = denoising_parameters.n_ext;

	clock_t end = clock();
	PCDNumber denoising_time = PCDNumber(end - begin) / CLOCKS_PER_SEC;
	cout << "**************************************************" << endl;
	cout << "TIMING" << endl;
	cout << "**************************************************" << endl;
	cout << "Full denoising time: " << denoising_time << endl;
	cout << "--------------------------------------------------" << endl;
	cout << "Prepare time: " << "n_exec=" << ext_n_ops << " time=" << prepare_time << " (" << (prepare_time / denoising_time)*100.0 << "%)" << endl;
	cout << "An. Neigh. time: " << "n_exec=" << ext_n_ops << " time=" << anisotropic_neighborhoods_time << " (" << (anisotropic_neighborhoods_time / denoising_time)*100.0 << "%)" << endl;
	cout << "Nor. filt. time: " << "n_exec=" << int_n_ops << " time=" << normal_filtering_time << " (" << (normal_filtering_time / denoising_time)*100.0 << "%)" << endl;
	cout << "Nor. corr. time: " << "n_exec=" << int_n_ops << " time=" << normal_correction_time << " (" << (normal_correction_time / denoising_time)*100.0 << "%)" << endl;
	cout << "Fe. class. time: " << "n_exec=" << int_n_ops << " time=" << feature_classification_time << " (" << (feature_classification_time / denoising_time)*100.0 << "%)" << endl;
	cout << "Po. update time: " << "n_exec=" << int_n_ops << " time=" << point_update_time << " (" << (point_update_time / denoising_time)*100.0 << "%)" << endl;
	cout << "**************************************************" << endl;

	return processed_point_cloud;
}
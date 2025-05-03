#ifndef UTIL_H
#define UTIL_H

#include "point_cloud.h"
#include "point_cloud_io.h"
#include "solver.h"
#include "convex_hull_2D.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <numeric>
#include <ilcplex/cplex.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/points/nn_search.h>
#include <geogram/basic/common.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/delaunay_2d.h>
#include <geogram/mesh/mesh.h>
#include <geogram/points/co3ne.h>
#include <geogram/mesh/mesh_geometry.h>
#include <stdio.h>
#include <stdlib.h> 

using namespace std;

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) 
{
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);
	stable_sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
	return idx;
}

PCDNumber compute_gaussian_weight(PCDNumber distance, PCDNumber sigma);

PCDNumber normal_distance(const PointCloud::Normal &n1, const PointCloud::Normal &n2);

PCDNumber point_to_line(PointCloud::Point & p, PointCloud::Point & p0, PointCloud::Point & v, PointCloud::Point & pb);

PCDNumber point_to_plane_distance(PointCloud::Point & p, PointCloud::Normal & plane_n, PointCloud::Point & plane_p);

PCDNumber point_to_plane_distance(PointCloud::Point & p, PointCloud::Normal & plane_n, PointCloud::Point & plane_p, PointCloud::Point & p_proj);

int plane_to_plane_intersection(PointCloud::Point & p1, PointCloud::Normal & n1, PointCloud::Point & p2, PointCloud::Normal & n2, PointCloud::Point & lp, PointCloud::Point & lv);

PCDNumber compute_rough_avg_edge_length(PointCloud & point_cloud);

void normalize_point_cloud(PointCloud & point_cloud, PCDNumber & scale_factor);

void retrieve_point_cloud(PointCloud & point_cloud, PCDNumber scale_factor);

void rescale_point_cloud(PointCloud & point_cloud, PCDNumber & scale_factor);

void compute_mesh_face_areas(PointCloud &mesh, vector<PCDNumber> &areas);

PCDNumber compute_mesh_vertex_area(PointCloud & mesh, PointCloud::VertexHandle vh, vector<PCDNumber> & areas);

pair<PCDNumber, PCDNumber> compute_2D_coordinate(PointCloud::Point & p, PointCloud::Point & plane_n, PointCloud::Point & plane_p, PointCloud::Point & plane_q);

PCDNumber point_to_plane_region_distance(PointCloud::Point & p, PointCloud::Point & plane_n, PointCloud::Point & plane_p, vector<PointCloud::Point> & considered_points);

PCDNumber angular_distance_on_tangent_plane(PointCloud::Point & normal, PointCloud::Point & feature_direction, PointCloud::Point & origin, PointCloud::Point & plane_normal, PointCloud::Point & query_point);

#endif // UTIL_H
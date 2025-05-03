#include "util.h"

PCDNumber compute_gaussian_weight(PCDNumber distance, PCDNumber sigma)
{
	return static_cast<PCDNumber>(exp(-0.5 * distance * distance / (sigma * sigma)));
}

PCDNumber normal_distance(const PointCloud::Normal &n1, const PointCloud::Normal &n2)
{
	return (n1 - n2).length();
}

PCDNumber point_to_line(PointCloud::Point & p, PointCloud::Point & p0, PointCloud::Point & v, PointCloud::Point & pb)
{
	PointCloud::Point w = p - p0;

	double c1 = w | v;
	double c2 = v | v;
	double b = c1 / c2;

	pb = p0 + b * v;
	return (p - pb).norm();
}

PCDNumber point_to_plane_distance(PointCloud::Point & p, PointCloud::Normal & plane_n, PointCloud::Point & plane_p)
{
	PCDNumber    sb, sn, sd;

	sn = -(plane_n | (p - plane_p));
	sd = (plane_n | plane_n);
	sb = sn / sd;

	PointCloud::Point b = p + sb * plane_n;
	return (p - b).length();
}

PCDNumber point_to_plane_distance(PointCloud::Point & p, PointCloud::Normal & plane_n, PointCloud::Point & plane_p, PointCloud::Point & p_proj)
{
	PCDNumber    sb, sn, sd;

	sn = -(plane_n | (p - plane_p));
	sd = (plane_n | plane_n);
	sb = sn / sd;

	PointCloud::Point b = p + sb * plane_n;
	p_proj = b;
	return (p - b).length();
}

int plane_to_plane_intersection(PointCloud::Point & p1, PointCloud::Normal & n1, PointCloud::Point & p2, PointCloud::Normal & n2, PointCloud::Point & lp, PointCloud::Point & lv)
{
	PointCloud::Point u = n1 % n2;
	PCDNumber ax = (u[0] >= 0 ? u[0] : -u[0]);
	PCDNumber ay = (u[1] >= 0 ? u[1] : -u[1]);
	PCDNumber az = (u[2] >= 0 ? u[2] : -u[2]);

	if ((ax + ay + az) < 0.00001) {
		PointCloud::Point v = p2 - p1;
		if ((n1 | v) == 0)
			return 1;
		else
			return 0;
	}
	int maxc;
	if (ax > ay)
	{
		if (ax > az)
			maxc = 1;
		else maxc = 3;
	}
	else
	{
		if (ay > az)
			maxc = 2;
		else maxc = 3;
	}

	PointCloud::Point iP;
	double d1, d2;
	d1 = -n1 | p1;
	d2 = -n2 | p2;

	switch (maxc)
	{
	case 1:
		iP[0] = 0;
		iP[1] = (d2*n1[2] - d1 * n2[2]) / u[0];
		iP[2] = (d1*n2[1] - d2 * n1[1]) / u[0];
		break;
	case 2:
		iP[0] = (d1*n2[2] - d2 * n1[2]) / u[1];
		iP[1] = 0;
		iP[2] = (d2*n1[0] - d1 * n2[0]) / u[1];
		break;
	case 3:
		iP[0] = (d2*n1[1] - d1 * n2[1]) / u[2];
		iP[1] = (d1*n2[0] - d2 * n1[0]) / u[2];
		iP[2] = 0;
	}

	lp = iP;
	lv = u;
	return 2;
}

PCDNumber compute_rough_avg_edge_length(PointCloud & point_cloud)
{
	int k = 10;

	GEO::initialize();
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("algo");

	double * vdata = new double[point_cloud.n_vertices() * 3];


	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		vdata[v_it->idx() * 3] = p[0];
		vdata[v_it->idx() * 3 + 1] = p[1];
		vdata[v_it->idx() * 3 + 2] = p[2];
	}

	GEO::NearestNeighborSearch_var NN = GEO::NearestNeighborSearch::create(3, "CNN");
	NN->set_points(point_cloud.n_vertices(), vdata);

	PCDNumber avg_edge_length = 0.0;
	PCDNumber avg_edge_length_count = 0.0;
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
	{
		PointCloud::Point p = point_cloud.point(*v_it);
		PCDNumber p_ptr[3]; p_ptr[0] = p[0]; p_ptr[1] = p[1]; p_ptr[2] = p[2];
		std::vector<GEO::index_t> neighs(k);
		std::vector<PCDNumber> sq_dists(k);
		NN->get_nearest_neighbors(k, p_ptr, neighs.data(), sq_dists.data());

		for (int i = 1; i < neighs.size(); i++)
		{
			avg_edge_length += sqrt(sq_dists[i]);
			avg_edge_length_count += 1.0;
		}
	}

	PCDNumber res = 0.1;
	if (avg_edge_length_count > 0)
		res = avg_edge_length / avg_edge_length_count;

	return res;
}

void normalize_point_cloud(PointCloud & point_cloud, PCDNumber & scale_factor)
{
	PCDNumber avg_edge_length = compute_rough_avg_edge_length(point_cloud);
	scale_factor = 1.0 / avg_edge_length;
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
		point_cloud.set_point(*v_it, (point_cloud.point(*v_it) * scale_factor));
}

void retrieve_point_cloud(PointCloud & point_cloud, PCDNumber scale_factor)
{
	PCDNumber temp = 1.0 / scale_factor;
	for (PointCloud::VertexIter v_it = point_cloud.vertices_begin(); v_it != point_cloud.vertices_end(); v_it++)
		point_cloud.set_point(*v_it, (point_cloud.point(*v_it) * temp));
}

void rescale_point_cloud(PointCloud & mesh, PCDNumber & scale_factor)
{
	for (PointCloud::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++)
		mesh.set_point(*v_it, (mesh.point(*v_it) * scale_factor));
}

void compute_mesh_face_areas(PointCloud &mesh, vector<PCDNumber> &areas)
{
	areas.resize(mesh.n_faces());

	for (PointCloud::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++)
	{
		vector<PointCloud::Point> point;
		point.resize(3); int index = 0;
		for (PointCloud::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++)
		{
			point[index] = mesh.point(*fv_it);
			index++;
		}
		PointCloud::Point edge1 = point[1] - point[0];
		PointCloud::Point edge2 = point[1] - point[2];
		PCDNumber S = 0.5f * (edge1 % edge2).length();
		areas[(*f_it).idx()] = S;
	}
}

PCDNumber compute_mesh_vertex_area(PointCloud & mesh, PointCloud::VertexHandle vh, vector<PCDNumber> & areas)
{
	PCDNumber area_sum = 0.0f;
	for (PointCloud::VertexFaceIter vf_iter = mesh.vf_iter(vh);vf_iter.is_valid();vf_iter++)
		area_sum += areas[vf_iter->idx()];
	return area_sum * (1.0f / 3.0f);
}

pair<PCDNumber, PCDNumber> compute_2D_coordinate(PointCloud::Point & p, PointCloud::Point & plane_n, PointCloud::Point & plane_p, PointCloud::Point & plane_q)
{
	PointCloud::Point neigh_point = p;
	PointCloud::Point vec = neigh_point - plane_p;
	PointCloud::Point proj_neigh_point;
	point_to_plane_distance(neigh_point, plane_n, plane_p, proj_neigh_point);
	PointCloud::Point centered_proj_neigh_point = proj_neigh_point - plane_p;
	PointCloud::Point u = (plane_q - plane_p).normalize_cond();
	PointCloud::Point v = (plane_n % u).normalize_cond();
	double u_neigh = centered_proj_neigh_point | u;
	double v_neigh = centered_proj_neigh_point | v;
	return pair<PCDNumber, PCDNumber>(u_neigh, v_neigh);
}

PCDNumber point_to_plane_region_distance(PointCloud::Point & p, PointCloud::Point & plane_n, PointCloud::Point & plane_p, vector<PointCloud::Point> & considered_points)
{
	if (considered_points.size() < 1)
		return 0;
	else if (considered_points.size() == 1)
		return (p - considered_points[0]).norm();
	else
	{
		Point2DCHCloud ch_pcl;

		PointCloud::Point rnp;
		point_to_plane_distance(considered_points[1], plane_n, plane_p, rnp);
		for (int i = 0; i < considered_points.size(); i++)
		{
			pair<PCDNumber, PCDNumber> coord_2D = compute_2D_coordinate(considered_points[i], plane_n, plane_p, rnp);
			ch_pcl.push_back(Point2DCH(coord_2D.first, coord_2D.second));
		}
		vector<int> ch_indices = convexHull(ch_pcl);

		pair<PCDNumber, PCDNumber> p_coord_2D = compute_2D_coordinate(p, plane_n, plane_p, rnp);

		Point2DCH p_proj(p_coord_2D.first, p_coord_2D.second);
		double min_d = MAX_PCDNUMBER;
		double cx = 0, cy = 0;
		bool lies_on_polygon = true;
		for (int i = 0; i < ch_indices.size(); i++)
		{
			Point2DCH p0 = ch_pcl[ch_indices[i]];
			Point2DCH p1 = ch_pcl[ch_indices[(i + 1) % ch_indices.size()]];
			cx += p0.x;
			cy += p0.y;
			double d = CH_point_to_line_distance(p_proj, p0, p1);
			if (d < min_d)
				min_d = d;

			int side = getSide(p0, p1, p_proj);
			if (!(side > 0))
				lies_on_polygon = false;
		}
		cx /= double(ch_indices.size());
		cy /= double(ch_indices.size());
		Point2DCH centroid(cx, cy);
		double cent_d = (p_proj - centroid).norm();
		if (lies_on_polygon)
			return 0.0;
		else return min_d;
		//return 0;
	}
}

PCDNumber angular_distance_on_tangent_plane(PointCloud::Point & normal, PointCloud::Point & feature_direction, PointCloud::Point & origin, PointCloud::Point & plane_normal, PointCloud::Point & query_point)
{
	PointCloud::Point vec = normal % feature_direction;
	vec = vec.normalize_cond();
	PointCloud::Point query_point_proj;
	PointCloud::Point vec_proj;
	point_to_plane_distance(query_point, plane_normal, origin, query_point_proj);
	PointCloud::Point query_vec = (query_point_proj - origin).normalize_cond();
	point_to_plane_distance(origin + vec, plane_normal, origin, vec_proj);
	vec_proj = vec_proj - origin;
	vec_proj = vec_proj.normalize_cond();
	PCDNumber angular_distance = abs(acos(abs(vec_proj | query_vec)));
	return angular_distance;
}
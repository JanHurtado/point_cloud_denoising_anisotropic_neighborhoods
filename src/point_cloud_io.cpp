#include "point_cloud_io.h"

bool import_point_cloud(PointCloud & pc, string filename)
{
	OpenMesh::IO::Options opt;
	opt += OpenMesh::IO::Options::VertexNormal;
	if (!OpenMesh::IO::read_mesh(pc, filename, opt))
	{
		cout << "IO error" << endl;
		return false;
	}
	return true;
}

bool export_point_cloud(PointCloud & pc, string filename)
{
	OpenMesh::IO::Options opt;
	opt += OpenMesh::IO::Options::VertexNormal;
	//opt += OpenMesh::IO::Options::VertexColor;
	return OpenMesh::IO::write_mesh(pc, filename, opt);
}
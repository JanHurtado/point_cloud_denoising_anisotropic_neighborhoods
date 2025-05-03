#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <iostream>
#include <queue>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/EigenValues>

using namespace std;

#define PI 3.14159265359
#define MAX_PCDNUMBER 999999.0f
#define MIN_PCDNUMBER -999999.0f

/** @addtogroup mesh_processing
  * @brief Mesh processing half-edge data structure definition (OpenMesh library).
  *
  * @{
  */

  /**
   * @brief The MyTraits struct - OpenMesh custom traits
   */
struct PointCloudTraits : OpenMesh::DefaultTraits
{
	// Let Point and Normal be a vector of AKNumbers
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;


	// The default 1D texture coordinate type is AKNumber.
	typedef double TexCoord1D;
	// The default 2D texture coordinate type is OpenMesh::Vec2f.
	typedef OpenMesh::Vec2f  TexCoord2D;
	// The default 3D texture coordinate type is OpenMesh::Vec3f.
	typedef OpenMesh::Vec3f  TexCoord3D;

	//enable standart properties
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);
	HalfedgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::PrevHalfedge);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color);
	EdgeAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Color);
};

/**
 * @brief AKMesh - triangular mesh definition
 */
typedef OpenMesh::TriMesh_ArrayKernelT<PointCloudTraits> PointCloud;

/**
 * @brief AKNumber - type definition for numbers used in mesh processing module
 */
typedef double PCDNumber;
typedef OpenMesh::Vec2d Point2D;

typedef Eigen::MatrixXf PCDMatrix;

/** @} */


#endif // POINT_CLOUD_H


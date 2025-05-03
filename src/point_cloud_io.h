#ifndef POINT_CLOUD_IO_H
#define POINT_CLOUD_IO_H

#include "point_cloud.h"

bool import_point_cloud(PointCloud & pc, string filename);

bool export_point_cloud(PointCloud & pc, string filename);

#endif // POINT_CLOUD_IO_H
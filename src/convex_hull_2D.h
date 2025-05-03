#ifndef CONVEX_HULL_2D_H
#define CONVEX_HULL_2D_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <list>
#include <set>

using namespace std;

const double CH_error = 0.00001;
const double CH_PI = 3.14159;
const double CH_max_value = 999999;
const double CH_min_value = -999999;

//point structure (used as vector too)
struct Point2DCH
{
	double x, y;
	Point2DCH() { x = 0; y = 0; }
	Point2DCH(double _x, double _y)
	{
		x = _x;
		y = _y;
	}
	Point2DCH& operator=(Point2DCH other)
	{
		x = other.x;
		y = other.y;
		return *this;
	}
	Point2DCH operator+(const Point2DCH & _p) const
	{
		Point2DCH res(x + _p.x, y + _p.y);
		return res;
	}
	Point2DCH operator-(const Point2DCH & _p) const
	{
		Point2DCH res(x - _p.x, y - _p.y);
		return res;
	}
	Point2DCH operator/(double scalar) const
	{
		Point2DCH res(x / scalar, y / scalar);
		return res;
	}
	Point2DCH operator*(double scalar) const
	{
		Point2DCH res(x * scalar, y * scalar);
		return res;
	}
	double operator*(const Point2DCH & _p) const
	{
		return x * _p.x + y * _p.y;
	}
	double norm()
	{
		return sqrt(x*x + y * y);
	}
};

//polygon centroid used to define counter-clockwise order
//extern Point2DCH CH_centroid;
//set of points
typedef vector<Point2DCH> Point2DCHCloud;

//defines the quadrant of a point
int quad(Point2DCH & p);

//sort by x-axis
bool sortX(pair<int, Point2DCH> & p1, pair<int, Point2DCH> & p2);

//sort by angle (counter-clockwise order)
class sortA
{
public:
	sortA(Point2DCH & _CH_centroid)
	{
		CH_centroid = _CH_centroid;
	}
	bool operator()(pair<int, Point2DCH> & p1, pair<int, Point2DCH> & p2)
	{
		Point2DCH _p1(p1.second.x - CH_centroid.x, p1.second.y - CH_centroid.y);
		Point2DCH _p2(p2.second.x - CH_centroid.x, p2.second.y - CH_centroid.y);

		int one = quad(_p1);
		int two = quad(_p2);

		if (one != two)
			return (one < two);
		return (_p1.y*_p2.x < _p2.y*_p1.x);
	}
private:
	Point2DCH CH_centroid;
};

//bool sortA(pair<int, Point2DCH> & p1, pair<int, Point2DCH> & p2);

//euclidean distance
double euclideanDistance(const Point2DCH & p1, const Point2DCH & p2);

//Given a rect defined by the points a and b, and a point p, this function returns the side where p is located (right, left or on the rect). Orientation.
int getSide(Point2DCH & a, Point2DCH & b, Point2DCH & p);

//extract a subset of vec where the initial index is a and the final index is b. a can be greater than b, assuming that the vector is a circular list.
vector<int> extractSubVector(vector<int> & vec, int a, int b);

//brute force Convex Hull. For all pair of points, we check if all the other points lie on one side with respect to the corresponding rect.
//returns the points in counter-clockwise order
vector<int> bruteHull(Point2DCHCloud & points, vector<int> & indices);

//merge function of Convex Hull
vector<int> convexHullMerge(Point2DCHCloud & _points, vector<int> & leftSide, vector<int> & rightSide);

//Convex Hull for a set of indices (Divide and Conquer)
vector<int> convexHull(Point2DCHCloud & _points, vector<int> indices);

//main convex hull function 
vector<int> convexHull(Point2DCHCloud & _points);

double CH_point_to_line_distance(Point2DCH & P, Point2DCH & P0, Point2DCH & P1);

#endif // CONVEX_HULL_2D_H
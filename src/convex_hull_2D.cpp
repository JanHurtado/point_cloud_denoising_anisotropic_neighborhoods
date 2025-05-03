#include "convex_hull_2D.h"

int quad(Point2DCH & p)
{
	if (p.x >= 0 && p.y >= 0)
		return 1;
	if (p.x <= 0 && p.y >= 0)
		return 2;
	if (p.x <= 0 && p.y <= 0)
		return 3;
	return 4;
}

//sort by x-axis
bool sortX(pair<int, Point2DCH> & p1, pair<int, Point2DCH> & p2) { return (p1.second.x < p2.second.x); }

//sort by angle (counter-clockwise order)


/*bool sortA(pair<int, Point2DCH> & p1, pair<int, Point2DCH> & p2)
{
	Point2DCH _p1(p1.second.x - CH_centroid.x, p1.second.y - CH_centroid.y);
	Point2DCH _p2(p2.second.x - CH_centroid.x, p2.second.y - CH_centroid.y);

	int one = quad(_p1);
	int two = quad(_p2);

	if (one != two)
		return (one < two);
	return (_p1.y*_p2.x < _p2.y*_p1.x);
}*/

//euclidean distance
double euclideanDistance(const Point2DCH & p1, const Point2DCH & p2)
{
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	return sqrt(dx*dx + dy * dy);
}

//Given a rect defined by the points a and b, and a point p, this function returns the side where p is located (right, left or on the rect). Orientation.
int getSide(Point2DCH & a, Point2DCH & b, Point2DCH & p)
{
	Point2DCH v = b - a;
	Point2DCH w = p - a;
	Point2DCH vOrt(-v.y, v.x);
	double val = vOrt * w;
	if (val == 0) return 0;
	else if (val > 0) return 1;
	else return -1;
}

//extract a subset of vec where the initial index is a and the final index is b. a can be greater than b, assuming that the vector is a circular list.
vector<int> extractSubVector(vector<int> & vec, int a, int b)
{
	vector<int> res;
	int n = vec.size();
	int i = a;
	for (int cont = 0; cont < n; cont++)
	{
		res.push_back(vec[i]);
		if (i == b) break;
		i = (i + 1) % n;
	}
	return res;
}

//brute force Convex Hull. For all pair of points, we check if all the other points lie on one side with respect to the corresponding rect.
//returns the points in counter-clockwise order
vector<int> bruteHull(Point2DCHCloud & points, vector<int> & indices)
{
	//compute the points that are part of the convex hull
	set<int> t_indices;
	for (int i = 0; i < indices.size(); i++)
	{
		for (int j = i + 1; j < indices.size(); j++)
		{
			int pos = 0, neg = 0;
			for (int k = 0; k < indices.size(); k++)
			{
				if (k != i && k != j)
				{
					Point2DCH a = points[indices[i]];
					Point2DCH b = points[indices[j]];
					Point2DCH p = points[indices[k]];
					if (getSide(a, b, p) >= 0)
						pos++;
					else
						neg++;
				}
			}
			if (pos == indices.size() - 2 || neg == indices.size() - 2)
			{
				t_indices.insert(indices[i]);
				t_indices.insert(indices[j]);
			}
		}
	}

	//sort the points in counter-clockwise order
	vector<pair<int, Point2DCH> >ret;
	for (set<int>::iterator it = t_indices.begin(); it != t_indices.end(); it++)
		ret.push_back(make_pair(*it, points[*it]));
	Point2DCH CH_centroid;
	CH_centroid.x = 0;
	CH_centroid.y = 0;
	int n = ret.size();
	for (int i = 0; i < n; i++)
	{
		CH_centroid.x += ret[i].second.x;
		CH_centroid.y += ret[i].second.y;
	}
	CH_centroid.x = CH_centroid.x / (double)n;
	CH_centroid.y = CH_centroid.y / (double)n;
	sort(ret.begin(), ret.end(), sortA(CH_centroid));
	vector<int> res;
	for (vector<pair<int, Point2DCH> >::iterator it = ret.begin(); it != ret.end(); it++)
	{
		res.push_back((*it).first);
	}
	return res;
}

//merge function of Convex Hull
vector<int> convexHullMerge(Point2DCHCloud & _points, vector<int> & leftSide, vector<int> & rightSide)
{
	int lowerTangentLeftSide;
	int lowerTangentRightSide;
	int upperTangentLeftSide;
	int upperTangentRightSide;

	//compute the corresponding most right and left values O(n)
	int mostRightLS = 0;
	double mostRightValueLS = CH_min_value;
	for (int i = 0; i < leftSide.size(); i++)
	{
		if (_points[leftSide[i]].x > mostRightValueLS)
		{
			mostRightValueLS = _points[leftSide[i]].x;
			mostRightLS = i;
		}
	}
	int mostLeftRS = 0;
	double mostLeftValueRS = CH_max_value;
	for (int i = 0; i < rightSide.size(); i++)
	{
		if (_points[rightSide[i]].x < mostLeftValueRS)
		{
			mostLeftValueRS = _points[rightSide[i]].x;
			mostLeftRS = i;
		}
	}

	//compute upper tangent O(n)
	bool flag = 0;
	int currentIdxLeftSide = mostRightLS;
	int currentIdxRightSide = mostLeftRS;
	int nLeftSide = leftSide.size();
	int nRightSide = rightSide.size();
	while (!flag)
	{
		flag = 1;
		while (getSide(_points[leftSide[currentIdxLeftSide]], _points[leftSide[(currentIdxLeftSide + 1) % nLeftSide]], _points[rightSide[currentIdxRightSide]]) < 0)
			currentIdxLeftSide = (currentIdxLeftSide + 1) % nLeftSide;

		while (getSide(_points[rightSide[currentIdxRightSide]], _points[rightSide[(nRightSide + currentIdxRightSide - 1) % nRightSide]], _points[leftSide[currentIdxLeftSide]]) > 0)
		{
			currentIdxRightSide = (nRightSide + currentIdxRightSide - 1) % nRightSide;
			flag = 0;
		}
	}
	upperTangentLeftSide = currentIdxLeftSide;
	upperTangentRightSide = currentIdxRightSide;

	//compute lower tangent O(n)
	currentIdxLeftSide = mostRightLS;
	currentIdxRightSide = mostLeftRS;
	flag = 0;
	while (!flag)
	{
		flag = 1;
		while (getSide(_points[rightSide[currentIdxRightSide]], _points[rightSide[(currentIdxRightSide + 1) % nRightSide]], _points[leftSide[currentIdxLeftSide]]) < 0)
			currentIdxRightSide = (currentIdxRightSide + 1) % nRightSide;

		while (getSide(_points[leftSide[currentIdxLeftSide]], _points[leftSide[(nLeftSide + currentIdxLeftSide - 1) % nLeftSide]], _points[rightSide[currentIdxRightSide]]) > 0)
		{
			currentIdxLeftSide = (nLeftSide + currentIdxLeftSide - 1) % nLeftSide;
			flag = 0;
		}
	}
	lowerTangentLeftSide = currentIdxLeftSide;
	lowerTangentRightSide = currentIdxRightSide;

	//combine the corresponding points preserving the counter-clockwise order
	vector<int> part1 = extractSubVector(leftSide, upperTangentLeftSide, lowerTangentLeftSide);
	vector<int> part2 = extractSubVector(rightSide, lowerTangentRightSide, upperTangentRightSide);
	part1.insert(part1.end(), part2.begin(), part2.end());

	return part1;
}

//Convex Hull for a set of indices (Divide and Conquer)
vector<int> convexHull(Point2DCHCloud & _points, vector<int> indices)
{
	if (indices.size() <= 6)
	{
		return bruteHull(_points, indices);
	}
	else
	{
		vector<int> leftIndices;
		for (int i = 0; i < indices.size() / 2; i++)
			leftIndices.push_back(indices[i]);
		vector<int> rightIndices;
		for (int i = indices.size() / 2; i < indices.size(); i++)
			rightIndices.push_back(indices[i]);
		vector<int> leftSide = convexHull(_points, leftIndices);
		vector<int> rightSide = convexHull(_points, rightIndices);
		return convexHullMerge(_points, leftSide, rightSide);
	}

}

//main convex hull function 
vector<int> convexHull(Point2DCHCloud & _points)
{
	//cout << "number of points: " << _points.size() << endl;
	vector<pair<int, Point2DCH> > data;
	for (int i = 0; i < _points.size(); i++)
		data.push_back(make_pair(i, _points[i]));
	//cout << "\t Sorting ..." << endl;
	sort(data.begin(), data.end(), sortX);
	Point2DCHCloud points;
	vector<int> indices;
	for (int i = 0; i < data.size(); i++)
	{
		indices.push_back(i);
		points.push_back(data[i].second);
		//cout<<i<<" "<<data[i].second.x<<" "<<data[i].second.y<<endl;
	}
	//cout << "\t Starting Convex Hull ..." << endl;
	//cout << "number of points: " << points.size() << endl;
	vector<int> convexHullLocalIndices = convexHull(points, indices);
	vector<int> convexHullGlobalIndices;
	for (int i = 0; i < convexHullLocalIndices.size(); i++)
		convexHullGlobalIndices.push_back(data[convexHullLocalIndices[i]].first);
	return convexHullGlobalIndices;
}

double CH_point_to_line_distance(Point2DCH & P, Point2DCH & P0, Point2DCH & P1)
{
	Point2DCH v = P1 - P0;
	Point2DCH w = P - P0;

	double c1 = w*v;
	double c2 = v*v;
	double b = c1 / (c2+0.0000001);

	Point2DCH Pb = P0 + v * b;
	return (P-Pb).norm();
}
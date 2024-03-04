#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include"point.h"
#include<vector>
#include<iostream>
#include"algorithm.h"
#include<set>
using namespace std;

/* Description: These functions are used to compute convex for a given point set
**
**                     Extreme edge, Gift wrapping, Quick hull, Gradam, Increment (by Bao)
**                     Divide and conquer (by https://www.tutorialspoint.com/convex-hull-using-divide-and-conquer-algorithm-in-cplusplus#:~:text=In%20this%20program%2C%20we%20will,to%20construct%20the%20convex%20hull.)
**                     
**					 1.ch_bykat()[4]; 
**                     2.ch_akl_toussaint()[5];
**                     3.ch_graham_andrew()[O(nlogn)  Graham-Andrew scan algorithm]
**                     4.ch_jarvis()[O(nh) Jarvis march algorithm]
**                     5.ch_eddy()[O(nh) Eddy's algorithm]
**                     6.ch_melkman()[linear time Melkman's algorithm] (this function is designed for handling polygon input formats)
**                     (by CGAL)
**                     
**Input          : Point set S
**Output       : Convex hull CH(S)
**History      : 02/26/24 Created by Xiyu Bao [bxy1310@gmail.com]
**                                    Add Extreme edge, Gift wrapping, Quick hull
**                    04/03/24 Update by Xiyu Bao
**                                    Add Gradam, Increment, Divide and conquer and convex hull functions in CGAL
*/
//********************************************************************
void convex_hull_extreme_edge(std::vector<point> &S, std::vector<point> &CH) {

	if (S.size() < 3) return;

	std::vector<int> next;
	for (size_t i = 0; i < S.size(); ++i)
		next.push_back(-1);

	for (int i=0;i<S.size();++i)
	{
		for (int j = 0; j < S.size(); ++j) {
			if (i == j) continue;
			bool valid_extrem_edge = true;
			for (int k = 0; k < S.size(); ++k) {
				if (k == i || k == j) continue;
				if (!isPointLeftOfVector(S[k], S[i], S[j])) {
					valid_extrem_edge = false;
					break;
				}
			}
			if (valid_extrem_edge) {
				next[i] = j;
			}
		}
	}

	int startIndex = -1;
	for (int i = 0; i < next.size(); ++i) {
		if (next[i] != -1) {
			startIndex = i;
			break;
		}
	}
	
	int nowIndex = startIndex;
	do {
		if (next[nowIndex] == -1) {
			std::cout << "program error" << std::endl;
		}
		/*if (next[nowIndex] != startIndex)
			std::cout << nowIndex << "->";
		else
			std::cout << nowIndex << std::endl;*/
		CH.push_back(S[nowIndex]);
		nowIndex = next[nowIndex];
	} while (startIndex != nowIndex);
}

#if 1
inline int find_lowest_point(std::vector<point>& S) {
	int index = -1;
	float minY = MaxValue, maxX = -1;
	for (int i = 0; i < S.size(); ++i) {
		if (S[i].y < minY - myeps) {
			index = i;
			minY = S[i].y;
			maxX = S[i].x;
		}
		else if (S[i].y < minY + myeps) {
			if (S[i].x < maxX) {
				index = i;
				maxX = S[i].x;
			}
		}
	}
	return index;
}
#endif

void convex_hull_gift_wrapping(std::vector<point>& S, std::vector<point>& CH) {
	if (S.size() < 3) return;
	int lowestIndex = find_lowest_point(S);
	int i = lowestIndex;
	point startEdge(1, 0);

	do {
		double theta_i = M_PI;
		int k;
		
		for (int j = 0; j < S.size(); ++j) {
			if (i == j) continue;
			double theta_j = counterClockwiseAngle(startEdge, S[j] - S[i]);
			/*if (S[i].x == 413 && S[i].y == 0)
				std::cout << S[j].x << "," << S[j].y << "," << theta_j << std::endl;*/
			if (theta_i > theta_j) {
				theta_i = theta_j;
				k = j;
			}
		}
		CH.push_back(S[k]);
		//std::cout << k << std::endl;
		startEdge = S[k] - S[i];
		i = k;
	
	} while (i != lowestIndex);
}

#if 1
inline int find_left_down_point(std::vector<point>& S) {
	int index = -1;
	float maxY = -1, minX = MaxValue;
	for (int i = 0; i < S.size(); ++i) {
		if (S[i].x < minX - myeps) {
			index = i;
			maxY = S[i].y;
			minX = S[i].x;
		}
		else if (S[i].x <  minX + myeps) {
			if (S[i].y > maxY) {
				index = i;
				maxY = S[i].y;
			}
		}
	}
	return index;
}


inline int find_right_down_point(std::vector<point>& S) {
	int index = -1;
	float maxY = -1, maxX = -1;
	for (int i = 0; i < S.size(); ++i) {
		if (S[i].x > maxX + myeps) {
			index = i;
			maxY = S[i].y;
			maxX = S[i].x;
		}
		else if (S[i].x > maxX - myeps) {
			if (S[i].y > maxY) {
				index = i;
				maxY = S[i].y;
			}
		}
	}
	return index;
}

inline int find_left_top_point(std::vector<point>& S) {
	int index = -1;
	float minY = MaxValue, minX = MaxValue;
	for (int i = 0; i < S.size(); ++i) {
		if (S[i].x < minX - myeps) {
			index = i;
			minY = S[i].y;
			minX = S[i].x;
		}
		else if (S[i].x < minX + myeps) {
			if (S[i].y < minY) {
				index = i;
				minY = S[i].y;
			}
		}
	}
	return index;
}

inline void Cut(std::vector<point>& S, point x, point y, std::vector<point>& S1, std::vector<point>& S2)
{
	for (int i = 0; i < S.size(); ++i) {
		if (S[i].equals(x) || S[i].equals(y)) continue; 
		double angle = counterClockwiseAngle(y - x, S[i]-x);
		if (angle > M_PI + myeps && angle < 2 * M_PI - myeps) {
			S1.push_back(S[i]);
		}
		if (angle > myeps && angle < M_PI - myeps) {
			S2.push_back(S[i]);
		}
	}
}

inline void HullLoop(point a, point b, std::vector<point>&M, std::vector<point> &result)
{
	if (M.empty())
		return;
	int index = -1;
	double maxLength = -1;
	for (int i = 0; i < M.size(); ++i) {
		double dis = distance(a, b, M[i]);
		if (dis > maxLength) {
			maxLength = dis;
			index = i;
		}
	}
	point c = M[index];

	std::vector<point> A;
	std::vector<point> B;
	double right_angle;

	for (int i = 0; i < M.size(); ++i) {
		if (M[i].equals(a) || M[i].equals(b) || M[i].equals(c)) continue;
		right_angle = counterClockwiseAngle(c - a, M[i]-a);
		if (right_angle > M_PI + myeps && right_angle < 2 * M_PI - myeps) {
			A.push_back(M[i]);
		}
		right_angle = counterClockwiseAngle(b - c, M[i]-c);
		if (right_angle > M_PI + myeps && right_angle < 2 * M_PI - myeps) {
			B.push_back(M[i]);
		}
	}

	std::vector<point> R1, R2;
	HullLoop(a, c, A,R1);
	HullLoop(c, b, B, R2);

	result.insert(result.end(), R1.begin(), R1.end());
	result.push_back(c);
	result.insert(result.end(), R2.begin(), R2.end());
}
#endif

void convex_hull_quick_hull(std::vector<point>& S, std::vector<point>& CH) {
	if (S.size() < 3) return;
	int x = find_right_down_point(S);
	int y = find_left_top_point(S);
	std::vector<point> S1;
	std::vector<point> S2;
	Cut(S, S[x], S[y], S1, S2);
	std::vector<point> R1, R2;
	CH.push_back(S[x]);
	HullLoop(S[x], S[y], S1, R1);
	CH.insert(CH.end(), R1.begin(), R1.end());
	CH.push_back(S[y]);
	HullLoop(S[y], S[x], S2, R2);
	CH.insert(CH.end(), R2.begin(), R2.end());
	//for (int i = 0; i < CH.size(); ++i) {
	//	std::cout << CH[i].x << "," << CH[i].y << ",";
	//}std::cout << std::endl;
}

#if 1
class indexTheta
{
public:
	indexTheta(int index, double theta) {
		this->index = index;
		this->theta = theta;
	}
	int index;
	double theta;
	
	bool operator<(const indexTheta rh) const
	{
		// return id <= rh.id;
		return this->theta < rh.theta-myeps;
	}
};
#endif

#include<stack>
#include<algorithm>
void convex_hull_improved_graham(std::vector<point>& S, std::vector<point>& CH)
{
	if (S.size() < 3)
		return;
	int p = find_right_down_point(S);
	std::vector<indexTheta> vertices;
	for (int i = 0; i < S.size(); ++i) {
		if (i == p) continue;
		double theta = counterClockwiseAngle(point(1, 0), S[i] - S[p]);
		vertices.push_back(indexTheta(i, theta));
	}

	std::sort(vertices.begin(), vertices.end());
	for (int i = 0; i < vertices.size()-1; ++i) {
		if (abs(vertices[i].theta - vertices[i + 1].theta) < myeps) {
			if (distance(S[p], S[vertices[i].index]) < distance(S[p], S[vertices[i + 1].index]))
			{
				vertices.erase(vertices.begin() + i);
				i--;
			}
			else {
				vertices.erase(vertices.begin() + i + 1);
			}
		}
	}

	std::vector<int> stack;
	stack.push_back(p);
	stack.push_back(vertices[0].index);
	stack.push_back(vertices[1].index);
	int top = 2;
	int i = 2;
	
	while (i < vertices.size()) {
		if (top == 0) {
			stack.push_back(vertices[i].index);
			top++;
			i++;
		}
		double theta = counterClockwiseAngle(
			S[stack[top]] - S[stack[top-1]],
			S[vertices[i].index] - S[stack[top - 1]]);
		if (theta > myeps&&theta<M_PI-myeps) {
			stack.push_back(vertices[i].index);
			top++;
			i++;
		}
		else {
			stack.erase(stack.begin() + top);
			top--;
		}
	}

	for (int i = 0; i <= top; ++i) {
		CH.push_back(S[stack[i]]);
	}
}

#if 1

void tangentPoint(std::vector<point>& CH, point pk, int& left, int& right)
{
	for (int i = 0; i < CH.size(); ++i)
	{
		point last = CH[(CH.size() + i - 1) % CH.size()];
		point next = CH[(i + 1) % CH.size()];
		double angle1 = counterClockwiseAngle(pk - last, CH[i] - last);
		double angle2 = counterClockwiseAngle(pk - CH[i], next - CH[i]);

		if (angle1 >= 0 && angle1 <= M_PI + myeps&&
			angle2>=M_PI-myeps&&angle2<2*M_PI)
			right = i;
		if (angle2 >= 0 && angle2 <= M_PI + myeps &&
			angle1 >= M_PI - myeps && angle1 < 2 * M_PI)
			left = i;
	}
}

void conv(std::vector<point>& CH, point pk)
{
	int left, right;
	tangentPoint(CH, pk, left, right);
	if (left < right) {
		for (int i = left+1; i < right; ++i) {
			CH.erase(CH.begin() + left+1);
		}
		CH.insert(CH.begin() + left + 1, pk);
	}
	if (left > right) {
		for (int i = left + 1; i < CH.size(); ++i) {
			CH.erase(CH.begin() + left+1);
		}
		for (int i = 0; i < right; ++i) {
			CH.erase(CH.begin());
		}
		CH.insert(CH.begin(), pk);
	}
}
#endif

void convex_hull_incremental(std::vector<point>& S, std::vector<point>& CH)
{
	if (S.size() < 3) return;

	//compute init convex hull
	for(int i=0;i<3;++i)
		CH.push_back(S[i]);
	double angle = counterClockwiseAngle(S[2] - S[0], S[1] - S[0]);
	if (angle <M_PI - myeps && angle >  myeps) {
		point temp = CH[1];
		CH[1] = CH[2];
		CH[2] = temp;
	}

	for (int i = 3; i < S.size(); ++i)
	{
		conv(CH, S[i]);
	}
}


#if 1
point mid;
// determines the quadrant of a point
// (used in compare())
int quad(point p)
{
	if (p.x >= 0 && p.y >= 0)
		return 1;
	if (p.x <= 0 && p.y >= 0)
		return 2;
	if (p.x <= 0 && p.y <= 0)
		return 3;
	return 4;
}


// Checks whether the line is crossing the polygon
int orientation(point a, point b,
	point c)
{
	int res = (b.y - a.y) * (c.x - b.x) -
		(c.y - b.y) * (b.x - a.x);

	if (res == 0)
		return 0;
	if (res > 0)
		return 1;
	return -1;
}

// compare function for sorting
bool compare(point p1, point q1)
{
	point p(p1.x - mid.x,
		p1.y - mid.y);
	point q(q1.x - mid.x,
		q1.y - mid.y);

	int one = quad(p);
	int two = quad(q);

	if (one != two)
		return (one < two);
	return (p.y * q.x < q.y * p.x);
}

// Finds upper tangent of two polygons 'a' and 'b'
// represented as two vectors.
vector<point> merger(vector<point > a,
	vector<point > b)
{
	// n1 -> number of points in polygon a
	// n2 -> number of points in polygon b
	int n1 = a.size(), n2 = b.size();

	int ia = 0, ib = 0;
	for (int i = 1; i < n1; i++)
		if (a[i].x > a[ia].x)
			ia = i;

	// ib -> leftmost point of b
	for (int i = 1; i < n2; i++)
		if (b[i].x < b[ib].x)
			ib = i;

	// finding the upper tangent
	int inda = ia, indb = ib;
	bool done = 0;
	while (!done)
	{
		done = 1;
		while (orientation(b[indb], a[inda], a[(inda + 1) % n1]) >= 0)
			inda = (inda + 1) % n1;

		while (orientation(a[inda], b[indb], b[(n2 + indb - 1) % n2]) <= 0)
		{
			indb = (n2 + indb - 1) % n2;
			done = 0;
		}
	}

	int uppera = inda, upperb = indb;
	inda = ia, indb = ib;
	done = 0;
	int g = 0;
	while (!done)//finding the lower tangent
	{
		done = 1;
		while (orientation(a[inda], b[indb], b[(indb + 1) % n2]) >= 0)
			indb = (indb + 1) % n2;

		while (orientation(b[indb], a[inda], a[(n1 + inda - 1) % n1]) <= 0)
		{
			inda = (n1 + inda - 1) % n1;
			done = 0;
		}
	}

	int lowera = inda, lowerb = indb;
	vector<point> ret;

	//ret contains the convex hull after merging the two convex hulls
	//with the points sorted in anti-clockwise order
	int ind = uppera;
	ret.push_back(a[uppera]);
	while (ind != lowera)
	{
		ind = (ind + 1) % n1;
		ret.push_back(a[ind]);
	}

	ind = lowerb;
	ret.push_back(b[lowerb]);
	while (ind != upperb)
	{
		ind = (ind + 1) % n2;
		ret.push_back(b[ind]);
	}
	return ret;

}

// Brute force algorithm to find convex hull for a set
// of less than 6 points
vector<point> bruteHull(vector<point> a)
{
	// Take any pair of points from the set and check
	// whether it is the edge of the convex hull or not.
	// if all the remaining points are on the same side
	// of the line then the line is the edge of convex
	// hull otherwise not
	set<point >s;

	for (int i = 0; i < a.size(); i++)
	{
		for (int j = i + 1; j < a.size(); j++)
		{
			int x1 = a[i].x, x2 = a[j].x;
			int y1 = a[i].y, y2 = a[j].y;

			int a1 = y1 - y2;
			int b1 = x2 - x1;
			int c1 = x1 * y2 - y1 * x2;
			int pos = 0, neg = 0;
			for (int k = 0; k < a.size(); k++)
			{
				if (a1 * a[k].x + b1 * a[k].y + c1 <= 0)
					neg++;
				if (a1 * a[k].x + b1 * a[k].y + c1 >= 0)
					pos++;
			}
			if (pos == a.size() || neg == a.size())
			{
				s.insert(a[i]);
				s.insert(a[j]);
			}
		}
	}

	vector<point>ret;
	for (auto e : s)
		ret.push_back(e);

	// Sorting the points in the anti-clockwise order
	mid = { 0, 0 };
	int n = ret.size();
	for (int i = 0; i < n; i++)
	{
		mid.x += ret[i].x;
		mid.y += ret[i].y;
		ret[i].x *= n;
		ret[i].y *= n;
	}
	sort(ret.begin(), ret.end(), compare);
	for (int i = 0; i < n; i++)
		ret[i] = point(ret[i].x / n, ret[i].y / n);

	return ret;
}

// Returns the convex hull for the given set of points
vector<point> divide(vector<point> a)
{
	// If the number of points is less than 6 then the
	// function uses the brute algorithm to find the
	// convex hull
	if (a.size() <= 5)
		return bruteHull(a);

	// left contains the left half points
	// right contains the right half points
	vector<point>left, right;
	for (int i = 0; i < a.size() / 2; i++)
		left.push_back(a[i]);
	for (int i = a.size() / 2; i < a.size(); i++)
		right.push_back(a[i]);

	// convex hull for the left and right sets
	vector<point>left_hull = divide(left);
	vector<point>right_hull = divide(right);

	// merging the convex hulls
	return merger(left_hull, right_hull);
}
#endif

void convex_hull_divide_and_conquer(std::vector<point> S, std::vector<point>& CH)
{
	std::sort(S.begin(), S.end());
	CH = divide(S);
}


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/ch_graham_andrew.h>
#include<CGAL/ch_jarvis.h>
#include<CGAL/ch_eddy.h>
#include<CGAL/ch_melkman.h>
#include <vector>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::vector<Point_2> Points;

enum convex_hull_type { Default, Bradford, Bykat, Graham, Jarvis, Eddy, Melkman };
//Given a sequence of n input points with h extreme points, the function convex_hull_2() 
// uses either the output-sensitive O(nh) algorithm of Bykat [5] (a non-recursive version of the quickhull [4] algorithm) 
//or the algorithm of Akl and Toussaint, which requires O(nlogn) time in the worst case
//[4] C.Bradford Barber, David P.Dobkin, and Hannu Huhdanpaa.The Quickhull algorithm for convex hulls.ACM Trans.Math.Softw., 22(4) : 469–483, December 1996.
//[5] A.Bykat.Convex hull of a finite set of points in two dimensions.Inform.Process.Lett., 7 : 296–298, 1978.

//Other method :  1.ch_bykat()[4]; 
//                          2.ch_akl_toussaint()[5];
//                          3.ch_graham_andrew()[O(nlogn)  Graham-Andrew scan algorithm]
//                          4.ch_jarvis()[O(nh) Jarvis march algorithm]
//                          5.ch_eddy()[O(nh) Eddy's algorithm]
//                          6.ch_melkman()[linear time Melkman's algorithm]

//Warning: 6.ch_melkman() is used to compute the convex hull of **simple polygonal chains (or polygons)**
//                do not use it to directly compute the convex hull of point set.

void convex_hull_cgal(std::vector<point> S, std::vector<point>& CH, convex_hull_type convex_hull_algorithm=Default) {
	Points points, result;
	for (int i = 0; i < S.size(); ++i)
		points.push_back(Point_2(S[i].x, S[i].y));

	switch (convex_hull_algorithm)
	{
	case Default:
		CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(result)); break; //O(hn)/O(nlogn)
	case Bradford:
		CGAL::ch_bykat(points.begin(), points.end(), std::back_inserter(result)); break; //O(hn)
	case Bykat:
		CGAL::ch_akl_toussaint(points.begin(), points.end(), std::back_inserter(result)); break; //O(nlogn)
	case Graham:
		CGAL::ch_graham_andrew(points.begin(), points.end(), std::back_inserter(result)); break; //O(nlogn)
	case Jarvis:
		CGAL::ch_jarvis(points.begin(), points.end(), std::back_inserter(result)); break;//O(nh)
	case Eddy:
		CGAL::ch_eddy(points.begin(), points.end(), std::back_inserter(result)); break;//O(nh)
	case Melkman:
		CGAL::ch_melkman(points.begin(), points.end(), std::back_inserter(result)); break;//O(n)
	default:
		break;
	}

	for (int i = 0; i < result.size(); ++i)
		CH.push_back(point(result[i].x(), result[i].y()));	
}




#pragma region test

std::vector<point> generate_pointset_demo() {
	static std::vector<point> S{
		point{4, 2}, point{ 5,3 }, point{ 7,4 },
			point{ 8,8 }, point{ 6,7 }, point{ 3,9 }};
	return S;
}

std::vector<point> generate_pointset_large(int num) {
	static std::vector<point> S;
	S.clear();
	for (int i = 0; i < num; ++i) {
		int x = rand() % 1000;
		int y = rand() % 1000;
		S.push_back(point(x, y));
	}
	return S;
}

void test_convexhull(std::vector<point> S, bool CHvisible=true) {

	const int numberOfMethods = 13;
	std::vector<point> CHpool[numberOfMethods];
	std::vector<std::string> CHmethod =
	{ "-extreme edge","-gift wrapping", "-quick hull",
		"-improved graham", "-incremental","-divide and conquer",
		"-cgal default","-cgal bykat", "-cgal akl_toussaint", "-cgal graham_andrew",
	"-cgal jarvis", "-cgal eddy", "-cgal melkman" };
	time_t start[numberOfMethods], end[numberOfMethods];

	start[0] = clock();
	convex_hull_extreme_edge(S, CHpool[0]);
	end[0] = clock();
	start[1] = clock();
	convex_hull_gift_wrapping(S, CHpool[1]);
	end[1] = clock();
	start[2] = clock();
	convex_hull_quick_hull(S, CHpool[2]);
	end[2] = clock();
	start[3] = clock();
	convex_hull_improved_graham(S, CHpool[3]);
	end[3] = clock();
	start[4] = clock();
	convex_hull_incremental(S, CHpool[4]);
	end[4] = clock();
	start[5] = clock();
	if(S.size()<100)
		convex_hull_divide_and_conquer(S, CHpool[5]);
	end[5] = clock();
	for (int i = 0; i < numberOfMethods - 6; ++i)
	{
		start[i+6] = clock();
		convex_hull_cgal(S, CHpool[6 + i], convex_hull_type(i));
		end[i+6] = clock();
	}

	std::cout << "**************************************" << std::endl;
	for (int i = 0; i < numberOfMethods-1; ++i) {
		std::cout << "convex hull "<<CHmethod[i]<<" :" << CHpool[i].size() << std::endl;
		std::cout << "time cost :" << end[i] - start[i] <<"ms"<< std::endl;
		if (CHvisible)
		{
			int index = find_left_down_point(CHpool[i]);
			for (int j = index; j < CHpool[i].size(); ++j) {
				if (j != (CHpool[i].size()+index-1)%CHpool[i].size())
					std::cout << "(" << CHpool[i][j].x << "," << CHpool[i][j].y << "),";
				else
					std::cout << "(" << CHpool[i][j].x << "," << CHpool[i][j].y << ")" << std::endl;
			}
			for (int j = 0; j < index; ++j) {
				if (j != (CHpool[i].size() + index - 1) % CHpool[i].size())
					std::cout << "(" << CHpool[i][j].x << "," << CHpool[i][j].y << "),";
				else
					std::cout << "(" << CHpool[i][j].x << "," << CHpool[i][j].y << ")" << std::endl;
			}
		}
	}
	std::cout << "**************************************" << std::endl << std::endl;

}

#pragma region

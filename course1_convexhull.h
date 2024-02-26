#pragma once
#include"point.h"
#include<vector>
#include<iostream>
#include"algorithm.h"

/* Description: These three function are used to compute convex for a given point set
**                     The time complexity of them are O(n3), O(kn) and O(n), respectively
**Input          : Point set S
**Output       : Convex hull CH(S)
**History      : 02/26/24 Create by Xiyu Bao [bxy1310@gmail.com]
*/

void convex_hull_extreme_edge(std::vector<point> &S, std::vector<point> &CH) {

	if (S.empty()) return;

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

void convex_hull_gift_wrapping(std::vector<point>& S, std::vector<point>& CH) {
	int lowestIndex = find_lowest_point(S);
	int i = lowestIndex;
	point startEdge(1, 0);

	do {
		double theta_i = M_PI;
		int k;
		for (int j = 0; j < S.size(); ++j) {
			if (i == j) continue;
			double theta_j = counterClockwiseAngle(startEdge, S[j] - S[i]);
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

void convex_hull_quick_hull(std::vector<point>& S, std::vector<point>& CH) {
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

void test_convexhull() {
	std::vector<point> S{
		point{4, 2}, point{ 5,3 }, point{ 7,4 },
			point{ 8,8 }, point{ 6,7 }, point{ 3,9 }};
	std::vector<point> CH1;
	std::vector<point> CH2;
	std::vector<point> CH3;
	convex_hull_extreme_edge(S, CH1);
	convex_hull_gift_wrapping(S, CH2);
	convex_hull_quick_hull(S, CH3);

	std::cout << "convex hull extreme edge :" << CH1.size() << std::endl;
	for (int i = 0; i < CH1.size(); ++i) {
		if (i < CH1.size() - 1)
			std::cout << "(" << CH1[i].x << "," << CH1[i].y << "),";
		else
			std::cout << "(" << CH1[i].x << "," << CH1[i].y << ")" << std::endl;
	}
	std::cout << "convex hull gift wrappint :" << CH2.size() << std::endl;
	for (int i = 0; i < CH1.size(); ++i) {
		if (i < CH2.size() - 1)
			std::cout << "(" << CH2[i].x << "," << CH2[i].y << "),";
		else
			std::cout << "(" << CH2[i].x << "," << CH2[i].y << ")" << std::endl;
	}
	std::cout << "convex hull quick hull :"<< CH3.size() << std::endl;
	for (int i = 0; i < CH3.size(); ++i) {
		if (i < CH3.size() - 1)
			std::cout << "(" << CH3[i].x << "," << CH3[i].y << "),";
		else
			std::cout << "(" << CH3[i].x << "," << CH3[i].y <<")" << std::endl;
	}
}
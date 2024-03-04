#pragma once

//"algorithm.h" provides basic geometry functions

#include"point.h"
#include <cmath>

#define M_PI 3.1415926
#define myeps 1e-6
#define MaxValue 1000000

double counterClockwiseAngle(point a, point b);
double crossProduct(point p1, point p2);
double dotProduct(point a, point b);
bool isPointLeftOfVector(point p, point source, point target);

double counterClockwiseAngle(point a, point b)
{
    double dot = dotProduct(a, b);
    double cross = crossProduct(a, b);

    // 计算反三角函数来获取角度
    double angle = std::atan2(cross, dot);

    // 角度转换为正数
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    if (angle > 2 * M_PI - myeps)
        angle = 0;
    return angle;
}

double crossProduct(point p1, point p2) {
    return p1.x * p2.y - p1.y * p2.x;
}

double distance(point a, point b) {
    return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

double distance(point a, point b, point c)
{
    double x1 = a.x;
    double y1 = a.y;
    double x2 = b.x;
    double y2 = b.y;
    double x = c.x;
    double y = c.y;
    return abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - x1 * y2) / 
        sqrt((y2-y1)*(y2-y1)+(x2-x1)*(x2-x1));
}

// 函数用于计算两个向量point之间的点积
double dotProduct(point a, point b) {
    return a.x * b.x + a.y * b.y;
}

bool isPointLeftOfVector(point p, point source, point target)
{
    point vecA = p - source;
    point vecB = target - source;
    double angle = counterClockwiseAngle(vecB, vecA);
    return angle > -myeps && angle < M_PI + myeps;
}

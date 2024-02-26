#pragma once
#include<cmath>
class point
{
public:
	double x;
	double y;

    point(double x, double y) {
        this->x = x;
        this->y = y;
    }
    point operator+(const point& p) const {
        point q(this->x + p.x, this->y + p.y);
        return q;

    }
    point operator-(const point& p) const {
        point q(this->x - p.x, this->y - p.y);
        return q;
    }
    double operator*(const point& p) const {
        return  this->x * p.x + this->y * p.y;
    }
    bool equals(point p) {
        return abs(this->x - p.x) < 1e-6 && abs(this->y - p.y )< 1e-6;
    }
};

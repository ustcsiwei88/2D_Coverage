#ifndef VEC2H
#define VEC2H

#include<vector>
#include<cmath>

using namespace std;

struct vec2{
	double x, y;
	vec2(double x, double y):x(x),y(y){}
	vec2(){}
	vec2 operator * (const double & d) const{
		return vec2(x*d, y*d);
	}
	vec2 operator / (const double & d) const{
		return vec2(x/d, y/d);
	}
	vec2 operator - (const vec2 & ano) const{
		return vec2(x-ano.x, y-ano.y);
	}
	vec2 operator + (const vec2 & ano) const{
		return vec2(x+ano.x, y+ano.y);
	}
	bool operator == (const vec2 & ano){
		return x == ano.x && y==ano.y;
	}
};

typedef vector<vec2> vec2l;

inline double dis(vec2 p1, vec2 p2){
	return hypot(p1.x-p2.x, p1.y-p2.y);
}

inline double norm(vec2 p){
	return p.x*p.x + p.y*p.y;
}

struct polygon{
	vec2l pts;
	vector<double> dis; // edge distance
	double circumference = .0;
	void clear(){
		circumference = 0;
		dis.clear();
		pts.clear();
	}
	// polygon():circumference(0.0){}
};

#endif
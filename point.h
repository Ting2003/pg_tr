// ----------------------------------------------------------------//
// Filename : point.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Point class
// ----------------------------------------------------------------//
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added opeartor !=
//   * added this log

#ifndef __POINT_H__
#define __POINT_H__

#include <fstream>
using namespace std;

// coordinate
class Point{
public:
	Point();
	Point(long x, long y, long z);
	void set(long, long, long);
	friend bool operator == (Point & a, Point & b);
	friend bool operator != (Point & a, Point & b);
	friend ostream & operator << (ostream & , const Point & );
	long x,y,z;
};

inline void Point::set(long _x, long _y, long _z){
	x=_x;
	y=_y;
	z=_z;
}

#endif


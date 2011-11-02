// ----------------------------------------------------------------//
// Filename : map_mat.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// header file of MapMat class
// This class uses a map to store the element.
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:37:43 CST 2011
//   * created file

#ifndef __MAP_MAT_H__
#define __MAP_MAT_H__

#include <utility>
#include <iostream>
#include <map>
#include "vec.h"
using namespace std;

class Vec;
class MapMat{
	struct MapMatComp{
		bool operator() (const SizeTPair& l, const SizeTPair& r) const{
			if(l.first == r.first)
				return l.second < r.second;
			else
				return l.first < r.first;
		}
	};

	typedef map<SizeTPair, double, MapMatComp> PairMap;
public:
	MapMat();
	size_t size() const;
	friend ostream & operator <<(ostream & os, const MapMat & hm);
	double & operator [] (const SizeTPair & id);
	bool has_index(const SizeTPair & p) const;
	void diagonal_split(MapMat & L, MapMat & D, MapMat & U) const;

	// vector matrix multiplication
	friend Vec operator *(const Vec & x, const MapMat & A);
	friend Vec operator *(const MapMat & A, const Vec & x);

	// scale
	MapMat & operator *= (double scale);
	friend MapMat operator *(const MapMat & A, double scale);
	friend MapMat operator *(double scale, const MapMat & A);

	// addition
	MapMat & operator += (const MapMat & B);
	MapMat operator + (const MapMat & B) const;

	friend class Algebra;
private:
	PairMap m; // 
};

#endif

// ----------------------------------------------------------------//
// Filename : map_mat.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation file of MapMat class
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:37:43 CST 2011
//   * created file

#include <iostream>
#include <iomanip>
#include "global.h"
#include "map_mat.h"

MapMat::MapMat(){
	//PairMap::value_type v(SizeTPair(1,1), 10.0);
	//m.insert(v);
}

size_t MapMat::size() const{
	return m.size();
}

ostream & operator <<(ostream & os, const MapMat & hm){
	MapMat::PairMap::const_iterator it;
	for(it=hm.m.begin();it!=hm.m.end();++it){
		const SizeTPair & p = it->first;
		os<<scientific
		  <<"[" <<setw(OUTPUT_WIDTH_INDEX)<<p.first
		  <<", "<<setw(OUTPUT_WIDTH_INDEX)<<p.second
		  <<"] "<<setw(OUTPUT_WIDTH_FLOAT)<<it->second<<endl;
	}
	return os;
}

// Use a key (pair) to retrieve the value
MapMat::PairMap::mapped_type & MapMat::operator [] (const SizeTPair & id){
	return m[id];
}

bool MapMat::has_index(const SizeTPair & p) const{
	return m.find(p) != m.end();
}

// vector matrix multiplication
// x*A
Vec operator *(const Vec & x, const MapMat & A){
	size_t n = x.size();
	Vec b(n);
	MapMat::PairMap::const_iterator it;
	for(it=A.m.begin();it!=A.m.end();++it){
		const SizeTPair & p = it->first;
		size_t i = p.first;
		size_t j = p.second;
		double v = it->second;
		b.val[j] += v * x.val[i];
	}
	return b;
}

// A*x
Vec operator *(const MapMat & A, const Vec & x){
	size_t n = x.size();
	Vec b(n);
	MapMat::PairMap::const_iterator it;
	for(it=A.m.begin();it!=A.m.end();++it){
		const SizeTPair & p = it->first;
		size_t i = p.first;
		size_t j = p.second;
		double v = it->second;
		b.val[i] += v * x.val[j];
	}
	return b;
}

// scale
MapMat & MapMat::operator *= (double scale){
	MapMat::PairMap::iterator it;
	for(it=m.begin();it!=m.end();++it){
		it->second *= scale;
	}
	return (*this);
}

MapMat operator *(const MapMat & A, double scale){
	return MapMat(A) *= scale;
}

MapMat operator *(double scale, const MapMat & A){
	return MapMat(A) *= scale;
}

// addition
MapMat & MapMat::operator += (const MapMat & B){
	MapMat::PairMap::const_iterator it;
	for(it=B.m.begin();it!=B.m.end();++it){
		const SizeTPair & p = it->first;
		m[p] += it->second;
	}
	return (*this);
}

MapMat MapMat::operator + (const MapMat & B) const{
	return MapMat(*this) += B;
}

void MapMat::diagonal_split(MapMat & L, MapMat & D, MapMat & U) const{
       PairMap::const_iterator it;
       for(it=m.begin(); it!=m.end(); ++it){
               const SizeTPair & p = (*it).first;
               size_t i = p.first;
               size_t j = p.second;
               double x = (*it).second;  
               if(i == j)     D[p] += x;
               else if(i > j) L[p] += x;
               else           U[p] += x;
       }
}

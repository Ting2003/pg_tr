// ----------------------------------------------------------------//
// Filename : hash_mat.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation file of HashMat class
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:37:43 CST 2011
//   * created file

#include <iostream>
#include <iomanip>
#include "global.h"
#include "hash_mat.h"

HashMat::HashMat(){
	m.set_empty_key(SizeTPair(-1,-1)); // MUST used in google's version
	//PairMap::value_type v(SizeTPair(1,1), 10.0);
	//m.insert(v);
}

size_t HashMat::size() const{
	return m.size();
}

// note that the output is *unordered*
// For SGI hash_map, must use const_iterator
ostream & operator <<(ostream & os, const HashMat & hm){
	HashMat::PairMap::const_iterator it;
	for(it=hm.m.begin();it!=hm.m.end();++it){
		const HashMat::SizeTPair & p = it->first;
		os<<scientific
		  <<"[" <<setw(OUTPUT_WIDTH_INDEX)<<p.first
		  <<", "<<setw(OUTPUT_WIDTH_INDEX)<<p.second
		  <<"] "<<setw(OUTPUT_WIDTH_FLOAT)<<it->second<<endl;
	}
	return os;
}

// Use a key (pair) to retrieve the value
HashMat::PairMap::data_type & HashMat::operator [] (SizeTPair & id){
	return m[id];
}

bool HashMat::has_index(const SizeTPair & p) const{
	return m.find(p) != m.end();
}

// vector matrix multiplication
// x*A
Vec operator *(const Vec & x, const HashMat & A){
	size_t n = x.size();
	Vec b(n);
	HashMat::PairMap::const_iterator it;
	for(it=A.m.begin();it!=A.m.end();++it){
		const HashMat::SizeTPair & p = it->first;
		size_t i = p.first;
		size_t j = p.second;
		double v = it->second;
		b.val[j] += v * x.val[i];
	}
	return b;
}

// A*x
Vec operator *(const HashMat & A, const Vec & x){
	size_t n = x.size();
	Vec b(n);
	HashMat::PairMap::const_iterator it;
	for(it=A.m.begin();it!=A.m.end();++it){
		const HashMat::SizeTPair & p = it->first;
		size_t i = p.first;
		size_t j = p.second;
		double v = it->second;
		b.val[i] += v * x.val[j];
	}
	return b;
}

// scale
HashMat & HashMat::operator *= (double scale){
	HashMat::PairMap::iterator it;
	for(it=m.begin();it!=m.end();++it){
		it->second *= scale;
	}
	return (*this);
}

HashMat operator *(const HashMat & A, double scale){
	return HashMat(A) *= scale;
}

HashMat operator *(double scale, const HashMat & A){
	return HashMat(A) *= scale;
}

// addition
HashMat & HashMat::operator += (const HashMat & B){
	HashMat::PairMap::const_iterator it;
	for(it=B.m.begin();it!=B.m.end();++it){
		const HashMat::SizeTPair & p = it->first;
		m[p] += it->second;
	}
	return (*this);
}

HashMat HashMat::operator + (const HashMat & B) const{
	return HashMat(*this) += B;
}

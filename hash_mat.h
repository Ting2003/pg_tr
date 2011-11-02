// ----------------------------------------------------------------//
// Filename : hash_mat.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// header file of HashMat class
// This class uses a hash_map to store the element.
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:37:43 CST 2011
//   * created file


// TODO: can try boost::unordered_map and google-sparsehash
#ifndef __HASH_MAT_H__
#define __HASH_MAT_H__

#include <utility>
#include <iostream>
#include <ext/hash_map>
#include "google/dense_hash_map"
#include "vec.h"
using namespace std;
using namespace __gnu_cxx;
using google::dense_hash_map;

#define GOOGLE

// hash function class for pair<size_t, size_t>
// the hash function is copied from boost library
// reference:
// http://stackoverflow.com/questions/738054/hash-function-for-a-pair-of-long-long
// http://www.boost.org/doc/libs/1_37_0/doc/html/hash/combine.html
// http://www.boost.org/doc/libs/1_37_0/doc/html/hash/reference.html#boost.hash_combine
//
namespace __gnu_cxx{ 
template<> struct hash< pair<size_t,size_t> >{
		hash<size_t> sh;
		size_t operator()( const pair<size_t,size_t> & p) const {
			size_t seed = 0;
			seed ^= sh(p.first) + 0x9e3779b9 + (seed<<6) + (seed >> 2);
			seed ^= sh(p.second) + 0x9e3779b9 + (seed<<6) + (seed >> 2);
			return seed;
		}
	};
}//namespace

class Vec;
class HashMat{
public:
	typedef pair<size_t,size_t> SizeTPair;
	// choose the type of hash_map implementation here
#ifndef GOOGLE
	typedef hash_map<SizeTPair, double> PairMap;
#else
	typedef dense_hash_map< SizeTPair, double, hash<SizeTPair> > PairMap;
#endif

	HashMat();
	size_t size() const;
	friend ostream & operator <<(ostream & os, const HashMat & hm);
	PairMap::data_type & operator [] (pair<size_t,size_t> & id);
	bool has_index(const SizeTPair & p) const;

	// vector matrix multiplication
	friend Vec operator *(const Vec & x, const HashMat & A);
	friend Vec operator *(const HashMat & A, const Vec & x);

	// scale
	HashMat & operator *= (double scale);
	friend HashMat operator *(const HashMat & A, double scale);
	friend HashMat operator *(double scale, const HashMat & A);

	// addition
	HashMat & operator += (const HashMat & B);
	HashMat operator + (const HashMat & B) const;

	friend class Algebra;
private:
	//hash_map< pair<size_t,size_t>, double> m; // 
	PairMap m; // 
};

#endif

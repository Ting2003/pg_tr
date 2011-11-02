// ----------------------------------------------------------------//
// Filename : vec.h
// Author : Zigang Xiao <zxiao2@illinois.edu>
//
// header file of class Vec
// ----------------------------------------------------------------//
// - Zigang Xiao - Tue Jan 25 14:45:57 CST 2011
//   * created file
//
#ifndef __VEC_H__
#define __VEC_H__

#include <cstring>
#include <iostream>
#include "global.h"
#include "triplet.h"
using namespace std;

class Vec{
public:
	Vec(size_t size=0);
	Vec(const Vec & v);     // copy constructor
	// this is tricky: use the supplied *v as member 
	Vec(size_t size, double * v);
	~Vec();
	Vec & operator = (const Vec & v);
	void resize(size_t size);
	size_t size() const;
	double * get_val() const;
	void set_val(double * v);
	void fill_zeros();
	double & operator [] (const size_t idx);

	friend ostream & operator << (ostream & os, const Vec & vec);

	friend bool operator == (const Vec & l, const Vec & r);

	// matrix-vector multiplication
	friend Vec operator *(const Vec & x, const Matrix & A);
	friend Vec operator *(const Matrix & A, const Vec & x);

	// vector multiplication
	friend double operator * (const Vec & l, const Vec & r);

	// scale
	Vec & operator *= (double scale);
	friend Vec operator * (const Vec & vec, double scale);
	friend Vec operator * (double scale, const Vec & vec);

	// +
	Vec & operator += (const Vec & r);
	Vec operator + (const Vec & r) const;
	
	// -
	Vec & operator -= (const Vec & r);
	Vec operator - (const Vec & r) const;

	friend double distance_inf(Vec & v, Vec & w);

	friend class Algebra;
private:
	double * val;
	size_t n;
};

#endif

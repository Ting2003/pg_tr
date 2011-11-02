// ----------------------------------------------------------------//
// Filename : vec.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//
// implentation file of class Vec
// ----------------------------------------------------------------//
// - Zigang Xiao - Tue Jan 25 14:45:57 CST 2011
//   * created file
//
#include <cstdio>
#include <cassert>
#include <iostream>
#include <iomanip>
#include "util.h"
#include "global.h"
#include "vec.h"
using namespace std;

// constructor
Vec::Vec(size_t size):n(size){
	//cout<<"constructor"<<endl;
	val = new double[n];
	memset(val, 0, sizeof(double) * n);
}

// private constructor: for run-time optimization
// for performing return value optimization
Vec::Vec(size_t size, double * v){
	//cout<<"private constructor"<<endl;
	n = size;
	val = new double[n];
	memcpy(val, v, sizeof(double) * n);
}

// copy constructor
Vec::Vec(const Vec & v):val(NULL),n(0){
	//cout<<"copy constructor"<<endl;
	(*this) = v;
}

// assignment operator
Vec & Vec::operator = (const Vec & v){
//	cout<<"assignment"<<endl;
	if( v.n > n ){// reallocate
		resize(v.n);
	}
	n = v.n;
	memcpy(val, v.val, sizeof(double)*n);

	return *this;
}

// desctructor
Vec::~Vec(){
//	cout<<"desctructor"<<endl;
	delete [] val;
}

void Vec::resize(size_t size){
//	cout<<"resize"<<endl;
	assert(size>0);
	delete [] val;
	val = new double[size];
	n = size;
}

// compare two Vec objects to see if they equal
bool operator == (const Vec & l, const Vec & r) {
	if(l.n == r.n){
		for(size_t i=0;i<l.n;i++)
			if( !feqn(l.val[i],r.val[i]) )  
				return false;
		return true;
	}
	else
		return false;
}

size_t Vec::size() const{return n;}

ostream & operator << (ostream & os, const Vec & vec){
	for(size_t i=0;i<vec.n;i++)
		os<<setw(OUTPUT_WIDTH_INDEX)<<i
		  <<setw(OUTPUT_WIDTH_FLOAT)<<vec.val[i]<<endl;
	return os;
}

double * Vec::get_val() const{return val;}

void Vec::set_val(double * v){
	delete [] val;
	val=v;
}

void Vec::fill_zeros(){
	memset(val, 0, sizeof(double)*n);
}

double & Vec::operator [] (const size_t idx){
	assert(idx<n);
	return val[idx];
}

///////////////////////////////////////////////////////////////////////////////

// compute the product of two vectors
double operator * (const Vec & l, const Vec & r){
	assert(l.size() == r.size());
	size_t size = l.size();
	double value=0;
	for(size_t i=0;i<size;i++) value += l.val[i]*r.val[i];
	return value;
}

// scale
Vec & Vec::operator *= (double scale){
	for(size_t i=0;i<size();i++)
		val[i]*=scale;
	return *this;

}

Vec operator * (const Vec & vec, double scale){
	return Vec(vec)*=scale;
}

Vec operator * (double scale, const Vec & vec){
	return Vec(vec)*=scale;
}

// +
Vec & Vec::operator += (const Vec & r){
	assert(this->size() == r.size());
	for(size_t i=0;i<this->size();i++)
		this->val[i] += r.val[i];
	return *this;

}

Vec Vec::operator + (const Vec & r) const{
	return Vec(*this) += r;
}

// -
Vec & Vec::operator -= (const Vec & r){
	assert(this->size() == r.size());
	for(size_t i=0;i<this->size();i++)
		this->val[i] -= r.val[i];
	return *this;

}

Vec Vec::operator - (const Vec & r) const{
	return Vec(*this) -= r;
}

double distance_inf(Vec & v, Vec & w){
	double max_diff = 0.0;
	assert(v.size() == w.size());
	for(size_t i=0;i<v.size();i++){
		double diff = fabs(v[i] - w[i]);
		if( max_diff < diff ) max_diff = diff;
	}
	return max_diff;
}

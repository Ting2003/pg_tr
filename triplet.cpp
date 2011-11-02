// ----------------------------------------------------------------//
// Filename : triplet.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//
// class definition of triplet
// ----------------------------------------------------------------//
// - Zigang Xiao - Fri Oct 29 16:43:17 CDT 2010
//   * created file
#include <iostream>
#include <iomanip>
#include <cassert>
#include "triplet.h"
#include "global.h"
#include "util.h"
using namespace std;

Triplet::Triplet(){
	//push_back(0,0,1.0);
}

ostream & operator <<(ostream & os, const Triplet & t){
	for(size_t i=0;i<t.size();i++){
		os<<scientific
		  <<"[" <<setw(OUTPUT_WIDTH_INDEX)<<t.Ti[i]
		  <<", "<<setw(OUTPUT_WIDTH_INDEX)<<t.Tj[i]
		  <<"] "<<setw(OUTPUT_WIDTH_FLOAT)<<t.Tx[i]<<endl;
	}
	return os;
}

// stupidly merge elements at the same position
void Triplet::merge(){
	for(size_t k=0;k<size();k++){
		for(size_t l=k+1;l<size();l++){
			if( Ti[k] == Ti[l] && Tj[k] == Tj[l] ){
				Tx[k] += Tx[l];
				Ti.erase(Ti.begin()+l);
				Tj.erase(Tj.begin()+l);
				Tx.erase(Tx.begin()+l);
			}
		}
	}
}

void Triplet::clear(){
	Ti.clear();
	Tj.clear();
	Tx.clear();
	// row = 0;
}
// insert a triplet 
// if index < 0, simply ignore
void Triplet::push_back(size_t i, size_t j, double x){
	Ti.push_back(i);
	Tj.push_back(j);
	Tx.push_back(x);
}

// return the number of elements in a triplet instance
size_t Triplet::size() const {return Ti.size();}

void Triplet::to_arrays(size_t * Ti, size_t * Tj, double * Tx) const{
	assert(Ti!=NULL && Tj!=NULL && Tx!=NULL);
	vector_to_array<size_t>(this->Ti, Ti);
	vector_to_array<size_t>(this->Tj, Tj);
	vector_to_array<double>(this->Tx, Tx);
}

// splits the matrix into three parts:
// D for diagonal, L for lower triangular, U for upper triangular
// Note that D, U don't include diagonal
void Triplet::diagonal_split(Triplet & L, Triplet & D, Triplet & U) const{
	size_t n=Tx.size();
	for(size_t k=0;k<n;k++){
		size_t i = Ti[k], j = Tj[k];
	        double x = Tx[k];
		if(i == j)     D.push_back(i,j,x);
		else if(i > j) L.push_back(i,j,x);
		else           U.push_back(i,j,x);
	}
}

////////////////////////////////////////////////////////////////////////////

// A: Sparse matrix in triplet format
// b: vector
// compute x*A=b, return b
Vec operator *(const Vec & x, const Triplet & A){
	size_t n = A.get_row();
	Vec b(n);
	for(size_t k=0;k<A.size();k++){
		size_t i=A.Ti[k];
		size_t j=A.Tj[k];
		double v=A.Tx[k];
		b.val[j] += v*x.val[i];
	}
	return b;
}

// A: Sparse matrix in triplet format
// b: vector
// compute A*x=b, return b
Vec operator *(const Triplet & A, const Vec & x){
	size_t n = A.get_row();
	Vec b(n);
	for(size_t k=0;k<A.size();k++){
		size_t i=A.Ti[k];
		size_t j=A.Tj[k];
		double v=A.Tx[k];
		b.val[i] += v * (x.val[j]);
	}
	return b;
}

// scale a matrix
Triplet & Triplet::operator *= (double scale){
	for(size_t k=0;k<Tx.size();k++)
		Tx[k] *= scale;
	return *this;
}

Triplet operator *(const Triplet & A, double scale){
	return Triplet(A) *= scale;
}

Triplet operator *(double scale, const Triplet & A){
	return Triplet(A) *= scale;
}

// addition
Triplet & Triplet::operator += (const Triplet & B){
	Ti.insert(Ti.end(), B.Ti.begin(), B.Ti.end());
	Tj.insert(Tj.end(), B.Tj.begin(), B.Tj.end());
	Tx.insert(Tx.end(), B.Tx.begin(), B.Tx.end());
	return *this;
}

Triplet Triplet::operator + (const Triplet & B) const{
	return Triplet(*this) += B;
}

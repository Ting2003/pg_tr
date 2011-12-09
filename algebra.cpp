// ----------------------------------------------------------------//
// Filename : algebra.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// linear algebra library
// ----------------------------------------------------------------//
// - Zigang Xiao - Wed Jan 26 18:15:02 CST 2011
//   * algebra::solve

#include <cassert>
#include "global.h"
#include "util.h"
#include "vec.h"
//#include "umfpack.h"
#include "algebra.h"

// deliver the address of x
void Algebra::solve_CK(Matrix & A, cholmod_factor *L, cholmod_dense *&x, cholmod_dense *b, 
  cholmod_common *cm){
	//cholmod_factor *L;
	cm->final_super = false;
	cm->final_asis = false;
	CK_decomp(A, L, cm);
	//cholmod_print_factor(L, "L", cm);
	//cholmod_print_common("CM", cm);
	// then solve
	//t1 = clock();
	x = cholmod_solve(CHOLMOD_A, L, b, cm);
	//t2 = clock();
	//clog<<"time for solving: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	//cholmod_print_dense(x, "x", cm);
	//cholmod_free_factor(&L, cm);
}
// doing cholesky decomposition
void Algebra::CK_decomp(Matrix &A, cholmod_factor *&L, cholmod_common *cm){
	// doing factorization first
	cholmod_triplet * T;
	size_t n_row = A.get_row();
	size_t n_col = A.get_row();
	size_t nnz = A.size();

	int *Ti;
	int *Tj;
	double *Tx;
	int stype = -1; // lower triangular storage
	//-1;// lower triangular storage
	T = cholmod_allocate_triplet(n_row, n_col, nnz, stype, 
			CHOLMOD_REAL, cm);
	Ti = static_cast<int *>(T->i);
	Tj = static_cast<int *>(T->j);
	Tx = static_cast<double *>(T->x);
	// copy data into T
	//if(nnz<THRESHOLD){
		for(size_t k=0;k<nnz;k++){
			Ti[k] = A.Ti[k];
			Tj[k] = A.Tj[k];
			Tx[k] = A.Tx[k];
		}
	//}
#if 0
	else{
		size_t k=0;
#pragma omp parallel for private(k)
		for(k=0;k<nnz;k++){
			Ti[k] = A.Ti[k];
			Tj[k] = A.Tj[k];
			Tx[k] = A.Tx[k];
		}	
	}
#endif
	T->nnz = nnz;
	//A.Ti.clear();
	//A.Tj.clear();
	//A.Tx.clear();
	cholmod_sparse * A_cholmod;
	A_cholmod = cholmod_triplet_to_sparse(T, nnz, cm);

        //cholmod_print_triplet(T, "T", cm);
	// free the triplet pointer
	cholmod_free_triplet(&T, cm);

	//cm->supernodal = -1;
	L = cholmod_analyze(A_cholmod, cm);
	//L->ordering = CHOLMOD_NATURAL;
	cholmod_factorize(A_cholmod, L, cm);
	//cholmod_print_factor(L, "L", cm);
	cholmod_free_sparse(&A_cholmod, cm);
}

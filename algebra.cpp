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
#include "umfpack.h"
#include "algebra.h"

// solve x for linear system Ax=b
// NOTE: UF_long and size_t must have the same size!
void Algebra::solve(const Matrix & A, const Vec & b, Vec & x){
	assert(x.size() == b.size());
	assert(sizeof(UF_long) == sizeof(size_t));
	clock_t t1,t2;

	size_t n = b.size();
	//Vec x(n);
	double * _x = x.val;
	double * _b = b.val;

	size_t n_row = n;
	size_t n_col = n;
	size_t nz = A.size();

	// NOTE: DO NOT MODIFY. size must be n_col+1, see UMFPACK manual
	UF_long * Ti = new UF_long[nz];
	UF_long * Tj = new UF_long[nz];
	double * Tx = new double[nz];
	A.to_arrays((size_t*)Ti,(size_t*)Tj,Tx);

	UF_long * Ap = new UF_long[n_col+1]; 
	UF_long * Ai = new UF_long[nz];
	double *Ax = new double [nz];

	int status;
	double Control [UMFPACK_CONTROL];
	umfpack_dl_defaults (Control) ;
	status = umfpack_dl_triplet_to_col(n_row, n_col, nz, Ti, Tj, Tx, 
			Ap, Ai, Ax, (UF_long *) NULL);

	if( status < 0 ) {
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_zi_triplet_to_col failed\n") ;
	}

	double *null = (double *) NULL;
	void *Symbolic, *Numeric;

	t1=clock();
	status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, 
			&Symbolic, Control, null); 
	t2=clock();
	//clog<<"Symbolic time = "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	if( status < 0 ){
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_symbolic failed\n") ;
	}

	t1=clock();
	status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, 
			&Numeric, Control, null) ;
	t2=clock();
	//clog<<"Numeric time = "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	if( status < 0 ){
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_numeric failed\n") ;
	}

	umfpack_dl_free_symbolic (&Symbolic) ;

	t1=clock();
	status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax, _x, _b, 
				Numeric, Control, null) ;
	t2=clock();
	//clog<<"Solve time = "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	if( status < 0 ){
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_solve failed\n") ;
	}
	umfpack_dl_free_numeric (&Numeric) ;

	delete [] Ti;
	delete [] Tj;
	delete [] Tx;
	delete [] Ax;
	delete [] Ai;
	delete [] Ap;

	//return x;
}

// Given a Matrix in triplet format, get its column compressed form representation
// Note that this function will allocate memory to Ap, Ai, Ax 
void Algebra::get_col_compressed(Matrix & A, size_t &count, UF_long *Ap, UF_long *Ai, double *Ax){
//	A.merge();
	size_t nz = A.size();
	size_t n_row = count;
	size_t n_col = count;

	// NOTE: DO NOT MODIFY. size must be n_col+1, see UMFPACK manual
	UF_long * Ti = new UF_long[nz];
	UF_long * Tj = new UF_long[nz];
	double * Tx = new double[nz];
	A.to_arrays((size_t*)Ti,(size_t*)Tj,Tx);

	int status;
	double Control [UMFPACK_CONTROL];
	umfpack_dl_defaults (Control) ;
	status = umfpack_dl_triplet_to_col(n_row, n_col, nz, Ti, Tj, Tx, 
			Ap, Ai, Ax, (UF_long *) NULL);
	if( status < 0 ){
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_triplet_to_col failed\n") ;
	}

	delete [] Ti;
	delete [] Tj;
	delete [] Tx;
}

// deliver the address of x
void Algebra::solve_CK(Matrix & A, cholmod_factor *L, cholmod_dense *&x, cholmod_dense *b, 
  cholmod_common *cm){
	//cholmod_factor *L;
	cm->final_super = false;
	cm->final_asis = false;
	double t1, t2;
	t1 = omp_get_wtime();
	CK_decomp(A, L, cm);
	t2 = omp_get_wtime();
	clog<<"time for decomp. "<<1.0*(t2-t1)<<endl;
	//cholmod_print_factor(L, "L", cm);
	//cholmod_print_common("CM", cm);
	// then solve
	t1 = omp_get_wtime();
	x = cholmod_solve(CHOLMOD_A, L, b, cm);
	t2 = omp_get_wtime();
	clog<<"time for solving: "<<1.0*(t2-t1)<<endl;
	//cholmod_print_dense(x, "x", cm);
	//cholmod_free_factor(&L, cm);
}

// Given column compressed form of matrix A
// perform LU decomposition and store the result in Numeric
// n is the dimension of matrix A
void Algebra::LU_decomposition(int n, UF_long * Ap, UF_long * Ai,
		double * Ax, void ** p_Numeric){
	int status;
	double Control [UMFPACK_CONTROL];
	umfpack_dl_defaults (Control) ;
	
	double *null = (double *) NULL;
	void * Symbolic;

	// perform ordering
	status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, 
			&Symbolic, Control, null); 
	if( status < 0 ){
		clog<<"lu decomp status: "<<status<<endl;
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_symbolic failed\n") ;
	}

	// LU decomposition
	status = umfpack_dl_numeric(Ap, Ai, Ax, Symbolic, 
			p_Numeric, Control, null) ;
	if( status < 0 ){
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_numeric failed\n") ;
	}

	umfpack_dl_free_symbolic (&Symbolic) ;
}

void Algebra::solve_numeric(UF_long *Ap, UF_long *Ai, double *Ax, Vec &x, Vec &b_new, void *Numeric){
	// the matrix needs to be initialized
	assert( Ap != NULL && Ai != NULL && Ax != NULL);

	int status;
	double Control [UMFPACK_CONTROL];
	umfpack_dl_defaults (Control) ;
	double *null = (double *) NULL;

	status = umfpack_dl_solve(UMFPACK_A, Ap, Ai, Ax,
		x.get_val(), b_new.get_val(), Numeric, 
		Control, null) ;
	if( status < 0 ){
		umfpack_dl_report_status (Control, status) ;
		report_exit("umfpack_dl_solve failed\n") ;
	}
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
	int stype = 1; // upper triangular storage
	//-1;// lower triangular storage
	T = cholmod_allocate_triplet(n_row, n_col, nnz, stype, 
			CHOLMOD_REAL, cm);
	Ti = static_cast<int *>(T->i);
	Tj = static_cast<int *>(T->j);
	Tx = static_cast<double *>(T->x);
	// copy data into T
	for(size_t k=0;k<nnz;k++){
		Ti[k] = A.Ti[k];
		Tj[k] = A.Tj[k];
		Tx[k] = A.Tx[k];
	}
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

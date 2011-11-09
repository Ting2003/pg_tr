// ----------------------------------------------------------------//
// Filename : circuit.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//          Ting Yu <tingyu1@illinois.edu>
//
// implementation file of circuit.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:35:56 CST 2011
//   * Add UMFPACK support
// - Zigang Xiao - Tue Jan 25 17:19:21 CST 2011
//   * added framework of PCG
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <utility>
#include <cassert>
#include <vector>
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"
using namespace std;

//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

vector<LAYER_DIR> Circuit::layer_dir(MAX_LAYER);

// constructor of Circuit class, name is optional
Circuit::Circuit(string _name):name(_name),
	circuit_type(UNKNOWN), VDD(0.0){
	// add ground node
	Node * gnd = new Node(string("0"), Point(-1,-1,-1));
	gnd->rep = gnd;
	this->add_node(gnd);

	for(int i=0;i<MAX_LAYER;i++)
		layer_dir[i]=NA;
	id_map = NULL;
	Lx = NULL;
	Li = NULL;
	Lp = NULL;
	Lnz = NULL;
}

// Trick: do not release memory to increase runtime
Circuit::~Circuit(){
	for(size_t i=0;i<nodelist.size();i++) delete nodelist[i];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		for(size_t j=0;j<ns.size();j++) delete ns[j];
	}
	free(id_map);
	Lx = NULL;
	Li = NULL;
	Lp = NULL;
	Lnz = NULL;
}

void Circuit::check_sys() const{
	clog<<"**** CHECKING SYSTEM ENVIRONMENT ****"<<endl;
	clog<<"* int size     = "<< sizeof(int)<<endl;
	clog<<"* long size    = "<< sizeof(long)<<endl;
	clog<<"* size_t size  = "<< sizeof(size_t)<<endl;
	clog<<"* UF_long size = "<< sizeof(UF_long)<<endl;
	clog<<"* Max nodelist = "<<(size_t)nodelist.max_size()<<endl;
	clog<<"****            END              ****"<<endl<<endl;
}

// functor to be used in STL sort function
// order: y > x > z > flag 
// input: two node a, b
// return true if a < b, false o/w
// note that ground note are put to last
bool compare_node_ptr(const Node * a, const Node * b){
	if( a->is_ground() ) return false;
	if (b->is_ground() ) return true;

	if( a->pt.y == b->pt.y ){
		if( a->pt.x == b->pt.x ){
			if( a->pt.z == b->pt.z ){
				return (a->isS() > b->isS());
			}
			else{
				return (a->pt.z > b->pt.z);// top down
			}
		}
		else
			return ( a->pt.x < b->pt.x );
	}
	else
		return (a->pt.y < b->pt.y);
}

void *Circuit::thread_task(){
   for(size_t i=start[my_id];i<end[my_id];i++){
      Node * p = nodelist[i];
      node_id[p] = i;
   }
   return 0;
}

// sort the nodes according to their coordinate 
void Circuit::sort_nodes(){
   sort(nodelist.begin(), nodelist.end(), compare_node_ptr);
   // update node id mapping, 
   // NOTE: ground node will be the last
   assign_task(NTHREADS, nodelist.size());
   start_thread();

   // update node id mapping, 
   // NOTE: ground node will be the last
#if 0
   for(i=0;i<nodelist.size();i++){
      Node * p = nodelist[i];
      node_id[p] = i;
   }
#endif 
}

string Circuit::get_name() const{return this->name;}

ostream & operator << (ostream & os, const NodePtrVector & nodelist){
	for(size_t i=0;i<nodelist.size();i++)
		os<<*nodelist[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const NetPtrVector & nets){
	for(size_t i=0;i<nets.size();i++)
		os<<*nets[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const Circuit & ckt){
	os<<"Circuit ["<<ckt.name<<"] info:"<<endl;

	os<<"==== Nodes ===="<<endl;
	os<<ckt.nodelist;

	os<<"==== Reps  ===="<<endl;
	os<<ckt.replist;

	return os;
}

void Circuit::print(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%s  %.5e\n", nodelist[i]->name.c_str(), 
				nodelist[i]->value);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Computation Functions

// initialization before solving the circuit
// 1. sort the nodes
// 2. set node representatives
// 3. find node in which block, update count
// 4. get representative lists
void Circuit::solve_init(){
        sort_nodes();
	size_t size = nodelist.size()-1; // exclude ground node!!
	Node *p=NULL;
	for(size_t i=0,nr=0;i<size;i++){
		p=nodelist[i];

		// set the representative
		Net * net = p->nbr[TOP];
		if( p->isS()== Y) VDD = p->get_value();

		// test short circuit
		if( p->isS() !=Y && // Y must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself

		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			//rep_id[p] = nr; // save the id
			p->rid = nr;
			++nr;
		}
		//clog<<"p->rep: "<<*p->rep<<endl;
	}
	//clog<<"nodelist.size: "<<nodelist.size()<<endl;
	//clog<<"replist.size: "<<replist.size()<<endl;
}

// count number of nodes that can be merged
void Circuit::count_merge_nodes(){
	size_t size = replist.size();
	size_t count = 0;
	for(size_t i=0;i<size;i++){
		Node * p = replist[i];
		if(p->nbr[TOP] ==NULL && p->nbr[BOTTOM] ==NULL){
			if((p->nbr[EAST] !=NULL && p->nbr[WEST] !=NULL) 
			&& (p->nbr[NORTH] ==NULL && p->nbr[SOUTH]==NULL))
				count++;
			else if((p->nbr[NORTH] !=NULL && p->nbr[SOUTH] !=NULL)			      &&(p->nbr[EAST]==NULL && p->nbr[WEST] ==NULL))
				count++;
		}
	}
	//clog<<"number of nodes can be merged is: "<<count<<endl;
}

void Circuit::solve(Tran &tran){
	solve_LU(tran);
}


// stamp the matrix and solve
void Circuit::solve_LU_core(Tran &tran){
   size_t n = replist.size();	// replist dosn't contain ground node
   if( n == 0 ) return;		// No node    
   cm = &c;
   cholmod_start(cm);
   cm->print = 5;
   b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
   bp = static_cast<double *> (b->x);

   Matrix A;
   stamp_by_set(A, bp);
   make_A_symmetric(bp);
   A.set_row(n);
   Algebra::solve_CK(A, L, x, b, cm);
   xp = static_cast<double *> (x->x);
   
   // print out dc solution
   //cout<<"dc solution. "<<endl;
   //cout<<nodelist<<endl;

   // A is already being cleared   
   if(nodelist.size()!=replist.size())
      assign_task(NTHREADS, n);
   start_ptr_memset();
   link_tr_nodes(tran);

   bnew = cholmod_zeros(n,1,CHOLMOD_REAL, cm);
   bnewp = static_cast<double *>(bnew->x);
 
   double time = 0;
   int iter = 0;
   stamp_by_set_tr(A, bp, tran);
   make_A_symmetric_tr(bp, xp, tran);
   
   stamp_current_tr(bp, time);
   
   clock_t t1, t2;
   t1 = clock();
   Algebra::CK_decomp(A, L, cm);
   Lp = static_cast<int *>(L->p);
   Lx = static_cast<double*> (L->x);
   Li = static_cast<int*>(L->i) ;
   Lnz = static_cast<int *>(L->nz);

   //cholmod_print_factor(L, "L", cm);
   
   t2 = clock();
   clog<<"decomp cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
   A.clear();
//#if 0 
   /*********** the following 2 parts can be implemented with pthreads ***/
   // build id_map immediately after transient factorization
   id_map = new int [n];
   cholmod_build_id_map(CHOLMOD_A, L, cm, id_map);

   temp = new double [n];
   // then substitute all the nodes rid
   start_ptr_assign_rid();
   start_ptr_assign_bp();

   start_ptr_assign_xp();
   start_ptr_assign_xp_b();

   delete [] temp;
   /*****************************************/ 
//#endif
   // bnewp[i] = bp[i]
   start_ptr_assign();

   //stamp_current_tr(bnewp, tran, time);
   modify_rhs_tr(bnewp, xp, tran);

   t1 = clock(); 
   solve_eq(xp);
   //x = cholmod_solve(CHOLMOD_A, L, bnew, cm);
   //xp = static_cast<double *>(x->x);
 
   t2 = clock();
   clog<<"time for solve_tr: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

   save_tr_nodes(tran, xp);
   time += tran.step_t;
   t1 = clock();
   // then start other iterations
   while(time < tran.tot_t){// && iter < 0){
       // bnewp[i] = bp[i];
       //for(size_t i=0;i<n;i++)
	//bnewp[i] = bp[i];
	start_ptr_assign();


      // only stamps if net current changes
      // set bp into last state
      //stamp_current_tr(bnewp, tran, time);
      stamp_current_tr_1(bp, bnewp, time);
     // get the new bnewp
      modify_rhs_tr(bnewp, xp, tran);
    
	
      solve_eq(xp);
      //x = cholmod_solve(CHOLMOD_A, L, bnew, cm); 
      //xp = static_cast<double *>(x->x);

      save_tr_nodes(tran, xp);
      time += tran.step_t;
      iter ++;
   }
   t2 = clock();
   clog<<"1000 iters cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

   release_tr_nodes(tran);
   cholmod_free_dense(&x, cm);
   cholmod_free_dense(&b, cm);
   cholmod_finish(&c);
}

// solve the node voltages using direct LU
void Circuit::solve_LU(Tran &tran){
	solve_init();
	solve_LU_core(tran);
}

// given vector x that obtained from LU, set the value to the corresponding
// node in nodelist
void Circuit::get_voltages_from_LU_sol(double * x){
   size_t i;
   for(i=0;i<nodelist.size()-1;i++){
      Node * node = nodelist[i];
      size_t id = node->rep->rid;  // get rep's id in Vec
      double v = x[id];		// get its rep's value
      node->value = v;

   }
}



// copy node voltages from the circuit to a Vec
// from = true then copy circuit to x
// else copy from x to circuit
void Circuit::copy_node_voltages(double * x, size_t &size, bool from){
	if( from == true ){
		for(size_t i=0;i<size;i++)
			x[i] = nodelist[i]->value;
	}
	else{
		for(size_t i=0;i<size;i++){
			Node * node = nodelist[i];
			size_t id = node_id[node->rep];
			double v = x[id];
			node->value = v;
		}
	}
}



// stamp the net in each set, 
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_by_set(Matrix & A, double * b){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			for(size_t i=0;i<ns.size();i++){
				assert( fzero(ns[i]->value) == false );
				stamp_resistor(A, ns[i]);
			}
			break;
		case CURRENT:
			for(size_t i=0;i<ns.size();i++)
				stamp_current(b, ns[i]);
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(A, b, ns[i]);
			}
			break;
		case CAPACITANCE:
			//for(size_t i=0;i<ns.size();i++)
				//stamp_capacitance_dc(A, ns[i]);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++){
				stamp_inductance_dc(A, b, ns[i]);	
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

// stamp transient current values into rhs
void Circuit::stamp_current_tr(double *b, double &time){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type == CURRENT)
			for(size_t i=0;i<ns.size();i++)
				stamp_current_tr_net(b, ns[i], time);
	}
}

// stamp transient current values into rhs
void Circuit::stamp_current_tr_1(double *bp, double *b, double &time){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type == CURRENT)
			for(size_t i=0;i<ns.size();i++)
				stamp_current_tr_net_1(bp, b, ns[i], time);
	}
}

// stamp the transient matrix
void Circuit::stamp_by_set_tr(Matrix & A, double *b, Tran &tran){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			for(size_t i=0;i<ns.size();i++){
				assert( fzero(ns[i]->value) == false );
				stamp_resistor_tr(A, ns[i]);
			}
			break;
		case CURRENT:
			//for(size_t i=0;i<ns.size();i++)
				//stamp_current_tr(b, ns[i]);
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD_tr(b, ns[i]);
			}
			break;
		case CAPACITANCE:
			for(size_t i=0;i<ns.size();i++)
				stamp_capacitance_tr(A, ns[i], tran);
			break;
		case INDUCTANCE:
			for(size_t i=0;i<ns.size();i++){
				stamp_inductance_tr(A, ns[i], tran);	
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

// update rhs by transient nets
void Circuit::modify_rhs_tr(double * b, double *x, Tran &tran){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_c_tr(ns[i], b, x, tran);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr(ns[i], b, x, tran);	
			}
		}
	}
}

void Circuit::stamp_resistor(Matrix & A, Net * net){
	//clog<<"net: "<<*net<<endl;
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
        if( !nk->is_ground()&& nk->isS()!=Y && 
          (nk->nbr[TOP]== NULL 
           || nk->nbr[TOP]->type != INDUCTANCE)) {
           A.push_back(k,k, G);

           //clog<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
           if(!nl->is_ground() &&(nl->nbr[TOP]==NULL || 
                 nl->nbr[TOP]->type != INDUCTANCE)&&(l > k)){
                    A.push_back(k,l,-G);

                  //clog<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
           }
        }

	if( !nl->is_ground() && nl->isS() !=Y && 
			(nl->nbr[TOP] ==NULL 
		||nl->nbr[TOP]->type != INDUCTANCE)) {
		A.push_back(l,l, G);
                //clog<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
		if(!nk->is_ground()&& (nk->nbr[TOP]==NULL ||
		  nk->nbr[TOP]->type != INDUCTANCE) && k > l){
			A.push_back(l,k,-G);
		//clog<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
                }
	}
}

// only stamp the resistor node connected to inductance
void Circuit::stamp_resistor_tr(Matrix & A, Net * net){
   double G;
   Node * nk = net->ab[0]->rep;
   Node * nl = net->ab[1]->rep;
   size_t k = nk->rid;
   size_t l = nl->rid;
   G = 1./net->value;

   if( nk->isS()!=Y && !nk->is_ground()&& 
     (nk->nbr[TOP]!=NULL && 
      nk->nbr[TOP]->type == INDUCTANCE)) {
      //clog<<"net: "<<*net<<endl;
      A.push_back(k,k, G);
      //clog<<"("<<k<<" "<<k<<" "<<G<<")"<<endl;
      if(!nl->is_ground() &&(nl->nbr[TOP]==NULL || 
           nl->nbr[TOP]->type != INDUCTANCE)){
         if(l > k){
            A.push_back(k,l,-G);
            //clog<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
         }
         else if(l < k){ 
            A.push_back(l, k, -G);
            //clog<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
         }
      }
   }

   if( nl->isS() !=Y && !nl->is_ground()&&
     (nl->nbr[TOP]!=NULL &&
      nl->nbr[TOP]->type == INDUCTANCE)) {

      //clog<<"net: "<<*net<<endl;
      A.push_back(l,l, G);
      //clog<<"("<<l<<" "<<l<<" "<<G<<")"<<endl;
      if(!nk->is_ground()&& (nk->nbr[TOP]==NULL ||
           nk->nbr[TOP]->type != INDUCTANCE)){
         if(k > l){
            A.push_back(l,k,-G);
            //clog<<"("<<l<<" "<<k<<" "<<-G<<")"<<endl;
         }
         else if(k < l){
            A.push_back(k, l, -G);
            //clog<<"("<<k<<" "<<l<<" "<<-G<<")"<<endl;
         }
      }
   }
}

void Circuit::stamp_inductance_dc(Matrix & A, double *b, Net * net){
	//clog<<"net: "<<*net<<endl;
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( nk->isS()!=Y && !nk->is_ground()) {
		A.push_back(k,k, 1);
		// general stamping
		if(!nl->is_ground())
		// A.push_back(k,l,-1);
		// make it symmetric
			b[k] = b[l];
		//clog<<"("<<k<<" "<<k<<" "<<1<<")"<<endl;
		//clog<<"("<<k<<" "<<l<<" "<<-1<<")"<<endl;
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		A.push_back(l,l, 1);
		if(!nk->is_ground())
		// general stamping
		// A.push_back(l,k,-1);
		b[l] = b[k];
		//clog<<"("<<l<<" "<<l<<" "<<1<<")"<<endl;
		//clog<<"("<<l<<" "<<k<<" "<<-1<<")"<<endl;
	}
}

// stamp dc cap, infinitive resistance
void Circuit::stamp_capacitance_dc(Matrix & A, Net * net){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( nk->isS()!=Y && !nk->is_ground()) {
		A.push_back(k,k, 0);
		if(!nl->is_ground())
			A.push_back(k,l,0);
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		A.push_back(l,l, 0);
		if(!nk->is_ground())
			A.push_back(l,k,0);
	}
}

// stamp inductance Geq = delta_t/(2L)
void Circuit::stamp_inductance_tr(Matrix & A, Net * net, Tran &tran){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = delta_t / (2*L)
	Geq = tran.step_t / (2*net->value);

	if( nk->isS()!=Y  && !nk->is_ground()) {
		// -1 is to clear formal inserted 1 at (k,k)
		A.push_back(k,k, Geq-1);
		//clog<<"("<<k<<" "<<k<<" "<<Geq-1<<")"<<endl;
		//clog<<nl->isS()<<endl;
		if(!nl->is_ground()&& nl->isS()!=Y && k<l){
			A.push_back(k,l,-Geq);
		        //clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		// -1 is to clear formal inserted 1 at (l,l)
		A.push_back(l,l, Geq-1);
		//clog<<"("<<l<<" "<<l<<" "<<Geq-1<<")"<<endl;
		if(!nk->is_ground() && nk->isS()!=Y && l<k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// stamp capacitance Geq = 2C/delta_t
void Circuit::stamp_capacitance_tr(Matrix &A, Net *net, Tran &tran){
	//clog<<"net: "<<*net<<endl;
	double Geq = 0;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	// Geq = 2*C / delta_t
	Geq = (2*net->value) / tran.step_t;
	//clog<<"C delta_t Geq: "<<net->value<<" "<<tran.step_t<<" "<<Geq<<endl;
	// Ieq = i(t) + 2*C / delta_t * v(t)

	if( nk->isS()!=Y  && !nk->is_ground()) {
		A.push_back(k,k, Geq);
		//clog<<"("<<k<<" "<<k<<" "<<Geq<<")"<<endl;
		if(!nl->is_ground()&& k < l){
			A.push_back(k,l,-Geq);
			//clog<<"("<<k<<" "<<l<<" "<<-Geq<<")"<<endl;
		}
	}

	if( nl->isS() !=Y && !nl->is_ground()) {
		A.push_back(l,l, Geq);
		//clog<<"("<<l<<" "<<l<<" "<<Geq<<")"<<endl;
		if(!nk->is_ground()&& l < k){
			A.push_back(l,k,-Geq);
			//clog<<"("<<l<<" "<<k<<" "<<-Geq<<")"<<endl;
		}
	}
}

// add Ieq into rhs
// Ieq = i(t) + 2*C / delta_t *v(t)
void Circuit::modify_rhs_c_tr(Net *net, double * rhs, double *x, Tran &tran){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	//clog<<"c net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
        // nk point to Z node
        if(nk->isS() != Z)
		swap<Node *>(nk, nl);
	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
	size_t k = nk->rid;
	size_t l = nl->rid;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) swap<Node *>(a, b);
	//clog<<"a, b: "<<*a<<" "<<*b<<endl;

	size_t id_a = a->rid;
	size_t id_b = b->rid;
	//i_t = (b->value - a->value) / r->value;
	i_t = (x[id_b] - x[id_a]) / r->value;
	//if(b->value != x[id_b] || a->value != x[id_a])
	   //cout<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;
	//clog<<"i_t: "<<i_t<<endl;
	//temp = 2*net->value / tran.step_t * 
		//(nk->value - nl->value);
       
        // push 2 nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
 #if 0
        if(iter ==0){
           pg.node_set_x.push_back(k);
           if(!nl->is_ground()) {
              //clog<<*nl<<" "<<l<<endl;
              pg.node_set_x.push_back(l);
           }
           else if(!b->is_ground()){
              //clog<<*b<<" "<<id_b<<endl;
              pg.node_set_x.push_back(id_b);
           }
        }
#endif
	if(nk->is_ground())
	 temp = 2*net->value / tran.step_t *
	       (0 - x[l]);
        else if(nl->is_ground()){
         temp = 2*net->value / tran.step_t *
	       (x[k] - 0);
        }
        else
         temp = 2*net->value / tran.step_t *
	       (x[k] - x[l]);
	//if(nk->value != x[k] || nl->value != x[l])
	   //cout<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;
	//clog<<"nk-nl "<<(nk->value - nl->value)<<" "<<2*net->value/tran.step_t<<" "<<temp<<endl;
	
	Ieq  = (i_t + temp);
	//clog<< "Ieq is: "<<Ieq<<endl;
	//clog<<"Geq is: "<<2*net->value / tran.step_t<<endl;
	if(!nk->is_ground()&& nk->isS()!=Y){
		 rhs[k] += Ieq;	// for VDD circuit
		//clog<<*nk<<" rhs +: "<<rhs[k]<<endl;
	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		 rhs[l] += -Ieq; 
		//clog<<*nl<<" rhs +: "<<rhs[l]<<endl;
	}
}
// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Circuit::modify_rhs_l_tr(Net *net, double *rhs, double *x, Tran &tran){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X) 
		swap<Node*>(nk, nl);
	size_t k = nk->rid;
	size_t l = nl->rid;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	//temp = tran.step_t / (2*net->value) * 
		//(nl->value - nk->value);
	temp = tran.step_t / (2*net->value) * 
		(x[l] - x[k]);
	//if(nk->value != x[k] || nl->value != x[l])
	   //clog<<"k, l, x_k, x_l: "<<nk->value<<" "<<nl->value<<" "<<
	     //x[k]<<" "<<x[l]<<endl;

	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	
	Net *r = nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	if(a->isS()!=X) swap<Node*>(a, b);
	size_t id_a = a->rid;
	size_t id_b = b->rid;
	i_t = (x[id_a] - x[id_b]) / r->value;
	//i_t = (a->value - b->value) / r->value;
        //if(b->value != x[id_b] || a->value != x[id_a])
	   //clog<<"a, b, x_a, x_b: "<<a->value<<" "<<b->value<<" "<<
	     //x[id_a]<<" "<<x[id_b]<<endl;

	//clog<<"resiste r: "<<*r<<endl;
	//clog<<*a<<" "<<*b<<endl;
	//clog<<"a, b, r, i_t: "<<a->value<<" "<<b->value<<" "<<
		//r->value<<" "<<i_t<<endl;
       
        // push inductance nodes into node_set_x
        //clog<<*nk<<" "<<k<<endl;
        //clog<<*b<<" "<<id_b<<endl;
 #if 0
        if(iter==0){
           pg.node_set_x.push_back(k);
           pg.node_set_x.push_back(id_b);
        }
#endif
	Ieq  = i_t + temp;
	//clog<<"Ieq: "<<Ieq<<endl;
	if(nk->isS() !=Y && !nk->is_ground()){
		 rhs[k] += Ieq; // VDD circuit
		//clog<<*nk<<" "<<rhs[k]<<endl;
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		 rhs[l] += -Ieq; // VDD circuit
		//clog<<*nl<<" "<<rhs[l]<<endl;
	}
}
// stamp a current source
void Circuit::stamp_current(double * b, Net * net){
	//clog<<"net: "<<*net<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;

	if( !nk->is_ground() && nk->isS()!=Y){// && 
		size_t k = nk->rid;
		b[k] += -net->value;
		//clog<<"b: "<<k<<" "<<-net->value<<endl;
	}
	if( !nl->is_ground() && nl->isS() !=Y){// &&
		size_t l = nl->rid;
		b[l] +=  net->value;
		//clog<<"b: "<<l<<" "<<-net->value<<endl;
	}
}

void Circuit::stamp_current_tr_net(double * b, Net * net, double &time){
	current_tr(net, time);
	//clog<<"net: "<<*net<<endl;
	//clog<<"current: "<<current<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground()&& nk->isS()!=Y) { 
		size_t k = nk->rid;
		//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
		b[k] += -net->value;//current;
		//clog<<"time, k, b: "<<time<<" "<<k<<" "<<b[k]<<endl;
	}
	if( !nl->is_ground() && nl->isS()!=Y) {
		size_t l = nl->rid;
		//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
		b[l] +=  net->value;// current;
		//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
	}
}

void Circuit::stamp_current_tr_net_1(double *bp, double * b, Net * net, double &time){
	double diff = 0;
	double current = net->value;
	current_tr(net, time);
	//if(time / 1e-11 == 22)
	//cout<<"time, co, cn: "<<time<<" "<<current<<" "<<net->value<<endl;
	// only stamps when net got a different current
	if(current != net->value){
		diff = net->value - current;
		
		//cout<<"time, old, new, diff: "<<time <<" "<<current<<" "<<net->value<<" "<<diff<<endl;
		//cout<<"net: "<<*net;
		//clog<<"current: "<<current<<endl;
		Node * nk = net->ab[0]->rep;
		Node * nl = net->ab[1]->rep;
		if( !nk->is_ground()&& nk->isS()!=Y) { 
			size_t k = nk->rid;
			//clog<<"node, rid: "<<*nk<<" "<<k<<endl;
			//clog<<"time, k, b bef: "<<time<<" "<<k<<" "<<b[k]<<endl;
			b[k] += -diff;//current;
			bp[k] = b[k];
			//clog<<"time, k, b: "<<time <<" "<<k<<" "<<b[k]<<endl;
		}
		if( !nl->is_ground() && nl->isS()!=Y) {
			size_t l = nl->rid;
			//clog<<"time, l, b bef: "<<time<<" "<<l<<" "<<b[l]<<endl;
			//clog<<"node, rid: "<<*nl<<" "<<l<<endl;
			b[l] +=  diff;// current;
			bp[l] = b[l];
			//clog<<"time, l, b: "<<time<<" "<<l<<" "<<b[l]<<endl;
		}
	}
}

// stamp a voltage source
void Circuit::stamp_VDD(Matrix & A, double * b, Net * net){
	// find the non-ground node
	//clog<<"net: "<<*net<<endl;
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	A.push_back(id, id, 1.0);
	//clog<<"push id, id, 1: "<<id<<" "<<id<<" "<<1<<endl;
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
		//clog<<"b: ="<<id<<" "<<net->value<<endl;
	}
	else{
		b[id] += net->value;
		//clog<<"b: +"<<id<<" "<<net->value<<endl;
	}
}

// stamp a voltage source
void Circuit::stamp_VDD_tr(double * b, Net * net){
	// find the non-ground node
	//clog<<"net: "<<*net<<endl;
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	//A.push_back(id, id, 1.0);
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
		//clog<<"b: ="<<id<<" "<<net->value<<endl;
	}
	else{
		b[id] += net->value;
		//clog<<"b: +"<<id<<" "<<net->value<<endl;
	}
}

// decide transient step current values
void Circuit::current_tr(Net *net, double &time){
	double slope = 0;
	double Tr = net->Tr;
	double PW = Tr + net->PW;
	double Tf = PW + net->Tf;
	double t_temp = time - net->TD;
	double t = fmod(t_temp, net->Period);
	if(time <= net->TD)
		net->value = net->V1;
	else if(t > 0 && t<= Tr){
		slope = (net->V2 - net->V1) / 
			(net->Tr);
		net->value = net->V1 + t*slope;
	}
	else if(t > Tr && t<= PW)
		net->value = net->V2;
	else if(t>PW && t<=Tf){
		slope = (net->V1-net->V2)/(net->Tf);
		net->value = net->V2 + slope*(t-PW);
	}
	else
		net->value = net->V1;
	//return current;
}

// assign value back to transient nodes
void Circuit:: save_tr_nodes(Tran &tran, double *x){
   size_t id=0;
   for(size_t j=0;j<tran.nodes.size();j++){
      if(tran.nodes[j].node != NULL ){
         id = tran.nodes[j].node->rep->rid;
         tran.nodes[j].value.push_back(x[id]);
      }
   }
}

// link transient nodes with nodelist
void Circuit:: link_tr_nodes(Tran &tran){
   for(size_t i=0;i<nodelist.size();i++){
      for(size_t j=0;j<tran.nodes.size();j++){
         if(tran.nodes[j].flag ==1) continue;
         if(nodelist[i]->name == tran.nodes[j].name){
            tran.nodes[j].node = nodelist[i];
            tran.nodes[j].flag =1;
            // record the id for tran.nodes
            break;
         }
      }
   }
}

void Circuit:: release_tr_nodes(Tran &tran){
   for(size_t j=0;j<tran.nodes.size();j++){
      if(tran.nodes[j].flag != -1){
         tran.nodes[j].node = NULL;
      }
   }
}




void Circuit::make_A_symmetric(double *b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p=NULL, *q=NULL, *r =NULL;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
           assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==X || (*it)->ab[1]->rep->isS()==X)) continue;
           // node p points to X node
           if((*it)->ab[0]->rep->isS()==X && ((*it)->ab[0]->rep->nbr[TOP]!=NULL && 
                (*it)->ab[0]->rep->nbr[TOP]->type==INDUCTANCE)){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==X && ((*it)->ab[1]->rep->nbr[TOP]!=NULL && 
                (*it)->ab[1]->rep->nbr[TOP]->type==INDUCTANCE)){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           
           r = p->nbr[TOP]->ab[0]->rep;
           if(r->isS()!=Y) 
              r = p->nbr[TOP]->ab[1]->rep;

           size_t id = q->rid;
           double G = 1.0 / (*it)->value;
           
           b[id] += r->value * G;
        }
}

void Circuit::make_A_symmetric_tr(double *b, double *x, Tran &tran){
	int type = INDUCTANCE;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
           if( (*it) == NULL ) continue;
           assert( fzero((*it)->value) == false );
           if(!((*it)->ab[0]->rep->isS()==Y || (*it)->ab[1]->rep->isS()==Y)) continue;
           //clog<<"net: "<<*(*it)<<endl;
           // node p points to Y node
           if((*it)->ab[0]->rep->isS()==Y){
              p = (*it)->ab[0]->rep; q = (*it)->ab[1]->rep;
           } 
           else if((*it)->ab[1]->rep->isS()==Y ){
              p = (*it)->ab[1]->rep; q = (*it)->ab[0]->rep;
           }           
           //clog<<"p and q: "<<*p<<" "<<*q<<endl;

           size_t id = q->rid;
           double G = tran.step_t / ((*it)->value*2);
           
           //b[id] += p->value * G;
           b[id] += x[p->rid] *G;
           //clog<<"stamp p->value, G, b: "<<p->value<<" "<<G<<" "<<b[id]<<endl;
        }
}

void Circuit::start_thread(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      pthread_create(&tid[i],NULL, Circuit::call_thread_task, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void Circuit::assign_task(int N_proc, int N){
   size_t num = N / N_proc;
   int res = N % N_proc;
   size_t  tot = 0;
   size_t num_temp = 0;
   for(int i=0;i<N_proc;i++){
      my_id =i;
      num_temp = num;
      if(i < res)
         num_temp++;
      start[i] = tot;
      tot += num_temp;
      end[i]= tot;
      //clog<<"proc, N, start, end: "<<i<<" "<<N<<" "<<start[i]<<" "<<end[i]<<endl;
   }
}

void* Circuit::call_thread_task(void *arg){
   return ((Circuit*)arg)->thread_task();
}

void Circuit::start_ptr_memset(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      my_id =i;
      pthread_create(&tid[i],NULL, Circuit::call_ptr_memset, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void *Circuit::ptr_memset(){
  for(size_t i=start[my_id];i<end[my_id];i++){
     bp[i]=0;
     //clog<<"i, my_id, bp: "<<i<<" "<<my_id<<" "<<bp[i]<<endl;
  }
  return 0;
}
void *Circuit::call_ptr_memset(void *arg){
   return ((Circuit*)arg)->ptr_memset();
}

void Circuit::start_ptr_assign(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      my_id =i;
      pthread_create(&tid[i],NULL, Circuit::call_ptr_assign, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void *Circuit::ptr_assign(){
  for(size_t i=start[my_id];i<end[my_id];i++){
     bnewp[i]=bp[i];
     //flag_col[i] = false;
     //clog<<"i, my_id, bp: "<<i<<" "<<my_id<<" "<<bp[i]<<endl;
  }
  return 0;
}
void *Circuit::call_ptr_assign(void *arg){
   return ((Circuit*)arg)->ptr_assign();
}

// xp[i] = temp[id_map[i]];
void Circuit::start_ptr_assign_xp_b(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      my_id =i;
      pthread_create(&tid[i],NULL, Circuit::call_ptr_assign_xp_b, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void *Circuit::ptr_assign_xp_b(){
  for(size_t i=start[my_id];i<end[my_id];i++){
	xp[i] = temp[id_map[i]];     
  }
  return 0;
}


void *Circuit::call_ptr_assign_xp_b(void *arg){
   return ((Circuit*)arg)->ptr_assign_xp_b();
}

// tempp[i] = xp[i];
void Circuit::start_ptr_assign_xp(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      my_id =i;
      pthread_create(&tid[i],NULL, Circuit::call_ptr_assign_xp, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void *Circuit::ptr_assign_xp(){
  for(size_t i=start[my_id];i<end[my_id];i++){
	temp[i] = xp[i];     
  }
  return 0;
}


void *Circuit::call_ptr_assign_xp(void *arg){
   return ((Circuit*)arg)->ptr_assign_xp();
}

// bp[i] = temp[id_map[i]];
void Circuit::start_ptr_assign_bp(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      my_id =i;
      pthread_create(&tid[i],NULL, Circuit::call_ptr_assign_bp, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void *Circuit::ptr_assign_bp(){
  for(size_t i=start[my_id];i<end[my_id];i++){
	bp[i] = temp[id_map[i]];     
  }
  return 0;
}


void *Circuit::call_ptr_assign_bp(void *arg){
   return ((Circuit*)arg)->ptr_assign_bp();
}

// reassign rid and copy bp into temp
void Circuit::start_ptr_assign_rid(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      my_id =i;
      pthread_create(&tid[i],NULL, Circuit::call_ptr_assign_rid, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void *Circuit::ptr_assign_rid(){
  for(size_t i=start[my_id];i<end[my_id];i++){
     int id = id_map[i];
     replist[id]->rid = i;
     temp[i] = bp[i];
     //clog<<"i, my_id, bp: "<<i<<" "<<my_id<<" "<<bp[i]<<endl;
  }
  return 0;
}


void *Circuit::call_ptr_assign_rid(void *arg){
   return ((Circuit*)arg)->ptr_assign_rid();
}

void Circuit::start_ptr_assign_1(){
//#if 0
   for(int i=0;i<NTHREADS;i++){
      my_id =i;
      pthread_create(&tid[i],NULL, Circuit::call_ptr_assign_1, this);
      pthread_join(tid[i], NULL);
   }
//#endif
}

void *Circuit::ptr_assign_1(){
  for(size_t i=start[my_id];i<end[my_id];i++){
     xp[i]=bnewp[i];
     //clog<<"i, my_id, bp: "<<i<<" "<<my_id<<" "<<bp[i]<<endl;
  }
  return 0;
}
void *Circuit::call_ptr_assign_1(void *arg){
   return ((Circuit*)arg)->ptr_assign_1();
}
/*void Circuit::parse_path_table(){
   // build up nodelist info
      Node_G *node;
      for(int i=0;i<replist.size();i++){
         node = new Node_G();
         node->value = i;
         pg.nodelist.push_back(node);
      }
}

void Circuit::build_path_graph(cholmod_factor *L){
   clock_t t1, t2;
   t1 = clock();
   build_FFS_path(L);
   t2 = clock();
   //clog<<"FFS path cost: "<<t2-t1<<endl;

   t1 = clock();
   build_FBS_path(L);
   t2 = clock();
   //clog<<"FBS path cost: "<<t2-t1<<endl;

   t1 = clock();
   // only keep the 2 paths, switch from List into array
   len_path_b = pg.path_FFS.get_size();
   len_path_x = pg.path_FBS.get_size();

   path_b = new int[len_path_b];
   path_x = new int [len_path_x];
   
   Node_G *nd;
   nd = pg.path_FFS.first;
   for(int i=0;i<len_path_b;i++){
      path_b[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FFS.destroy_list();

   nd = pg.path_FBS.first;
   for(int i=0;i<len_path_x;i++){
      path_x[i] = nd->value;
      if(nd->next != NULL)
         nd = nd->next;
   }
   pg.path_FBS.destroy_list();

   pg.nodelist.clear();
   pg.node_set_b.clear();
   pg.node_set_x.clear();
   t2 = clock();
   //clog<<"post cost: "<<t2-t1<<endl;
}

void Circuit::build_FFS_path(cholmod_factor *L){
   clock_t t1, t2;
   t1 = clock();
   parse_path_table(); 
   t2 = clock();
   //clog<<endl<<"parse path table cost: "<<t2-t1<<endl;

   t1 = clock();
   set_up_path_table(L);
   t2 = clock();
   //clog<<"set up path table cost: "<<t2-t1<<endl;

   t1 = clock();
   find_path(pg.node_set_b, pg.path_FFS);
   t2 = clock();
   //clog<<"find path cost: "<<t2-t1<<endl;

   t1 = clock();
   pg.path_FFS.assign_size();
   t2 = clock();
   //clog<<"assign size cost: "<<t2-t1<<endl;

   t1 = clock();
   for(int i=0;i<replist.size();i++)
      pg.nodelist[i]->flag = 0;
   t2 = clock();
   //clog<<"zeros flag cost: "<<t2-t1<<endl;
}

void Circuit::build_FBS_path(cholmod_factor *L){
  pg.nodelist.clear();
  parse_path_table();
  set_up_path_table(L);
  find_path(pg.node_set_x, pg.path_FBS);
   pg.path_FBS.assign_size();
}

void Circuit::set_up_path_table(cholmod_factor *L){
   size_t n = L->n;
   int *Lp, *Li, *Lnz;
   int p, lnz, s, e;
   Lp = static_cast<int *> (L->p);
   Lnz = (int *)L->nz;
   Li = static_cast<int *> (L->i);
   for(size_t i=0;i<n;i++){
      p = Lp[i];
      lnz = Lnz[i];

      s = Li[p];
      e = s;
      if(lnz >1) 
         e = Li[p+1];

      if(s<e)
         pg.nodelist[s]->next = pg.nodelist[e];
   }
}

void Circuit::find_path(vector<size_t> &node_set, List_G &path){
   Node_G* ne = pg.nodelist[pg.nodelist.size()-1];
   vector <Node_G *> insert_list;
   //clock_t s, e;
   //s = clock();
   sort(node_set.begin(), node_set.end());
   //for(int i=0;i<node_set.size();i++)
      //clog<<"i, node_set: "<<i<<" "<<node_set[i]<<endl;
   if(node_set.size() == 0) return;
   
   // build up the first path start with id = min 
   int id = node_set[0];
   do{
      path.add_node(pg.nodelist[id]);
      pg.nodelist[id]->flag =1;
      if(pg.nodelist[id]->next == NULL) break;
      pg.nodelist[id] = pg.nodelist[id]->next;
   }while(pg.nodelist[id]->value != ne->value);
   path.add_node(ne);
  // e = clock();
   //clog<<"first path cost: "<<1.0*(e-s)/CLOCKS_PER_SEC<<endl;

   //s = clock();
   for(size_t i=0; i<node_set.size();i++){
      int id = node_set[i];
      if(pg.nodelist[id]->flag == 1) continue;
      // stops at first place where flag = 0
      // clog<<"stops at i: "<<i<<endl;
      do{
         if(pg.nodelist[id]->flag ==0){
            insert_list.push_back(pg.nodelist[id]);
            pg.nodelist[id]->flag =1;
         }
         if(pg.nodelist[id]->next == NULL || 
           pg.nodelist[id]->next->flag ==1)
            break;
         // else clog<<"next node: "<<*pg.nodelist[id]->next;
         pg.nodelist[id] = pg.nodelist[id]->next;
      }while(pg.nodelist[id]->value != ne->value); 
   }
   //e = clock();
   //clog<<"store insert nodes cost: "<<1.0*(e-s)/CLOCKS_PER_SEC<<endl;

   //s = clock();
   //clog<<"insert_list.size: "<<insert_list.size()<<endl;
   sort(insert_list.begin(), insert_list.end(), compare_Node_G);
   //e = clock();
   //clog<<"sort nodes cost: "<<1.0*(e-s)/CLOCKS_PER_SEC<<endl;
   //for(int i=0;i<insert_list.size();i++)
      //clog<<"i, insert: "<<i<<" "<<*insert_list[i]<<endl;

   //clog<<"path: "<<&path<<endl;
   //s = clock();
   // p is the old pointer to the list
   // will be updated into new one
   Node_G *q, *p;
//#if 0
   for(size_t k=0;k<insert_list.size();k++){
      if(k ==0) p = path.first;
      else p = q;
      q = path.insert_node(insert_list[k], p);
   }
//#endif
   //e = clock();
   //clog<<"insert node cost: "<<1.0*(e-s)/CLOCKS_PER_SEC<<endl;

   //clog<<"path: "<<&path<<endl;
   //clog<<endl;
   insert_list.clear();
}
*/




void Circuit::solve_eq(double *X){
   int p, q, r, lnz, pend;
   int j, n = L->n ;
   // assign xp[i] = bnewp[i]
   for(int i=0;i<n;i++)
	X[i] = bnewp[i];
   //start_ptr_assign_1(); 

   //clock_t t1, t2;
   //t1 = clock();
   // FFS solve
   for (j = 0 ; j < n ; ){
      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
      {

         /* -------------------------------------------------------------- */
         /* solve with a single column of L */
         /* -------------------------------------------------------------- */

         double y = X [j] ;
         if(L->is_ll == true){
            X[j] /= Lx [p] ;
         }

         for (p++ ; p < pend ; p++)
         {
            X [Li [p]] -= Lx [p] * y ;
         }
         j++ ;	/* advance to next column of L */

      }
      //#if 0
      else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of two columns of L */
         /* -------------------------------------------------------------- */

         double y [2] ;
         q = Lp [j+1] ;
         if(L->is_ll == true){
            y [0] = X [j] / Lx [p] ;
            y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
            X [j  ] = y [0] ;
            X [j+1] = y [1] ;
         }

         else{
            y [0] = X [j] ;
            y [1] = X [j+1] - Lx [p+1] * y [0] ;
            X [j+1] = y [1] ;
         }
         for (p += 2, q++ ; p < pend ; p++, q++)
         {
            X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
         }
         j += 2 ;	    /* advance to next column of L */

      }
      else
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of three columns of L */
         /* -------------------------------------------------------------- */

         double y [3] ;
         q = Lp [j+1] ;
         r = Lp [j+2] ;
         if(L->is_ll == true){
            y [0] = X [j] / Lx [p] ;
            y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
            y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
            X [j  ] = y [0] ;
            X [j+1] = y [1] ;
            X [j+2] = y [2] ;
         }

         else{
            y [0] = X [j] ;
            y [1] = X [j+1] - Lx [p+1] * y [0] ;
            y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
            X [j+1] = y [1] ;
            X [j+2] = y [2] ;
         }
         for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
         {
            X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
         }
         j += 3 ;	    /* advance to next column of L */
      }
      //#endif
   }
   //t2 = clock();
   //clog<<"FFS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

   //t1 = clock();
   // FBS solve
   for(j = n-1; j >= 0; ){

      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      /* find a chain of supernodes (up to j, j-1, and j-2) */

      if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
      {

         /* -------------------------------------------------------------- */
         /* solve with a single column of L */
         /* -------------------------------------------------------------- */

         double d = Lx [p] ;
         if(L->is_ll == false)
            X[j] /= d ;
         for (p++ ; p < pend ; p++)
         {
            X[j] -= Lx [p] * X [Li [p]] ;
         }
         if(L->is_ll == true)
            X [j] /=  d ;
         j--;
      }
      else if (lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of two columns of L */
         /* -------------------------------------------------------------- */

         double y [2], t ;
         q = Lp [j-1] ;
         double d [2] ;
         d [0] = Lx [p] ;
         d [1] = Lx [q] ;
         t = Lx [q+1] ;
         if(L->is_ll == false){
            y [0] = X [j  ] / d [0] ;
            y [1] = X [j-1] / d [1] ;
         }
         else{
            y [0] = X [j  ] ;
            y [1] = X [j-1] ;
         }
         for (p++, q += 2 ; p < pend ; p++, q++)
         {
            int i = Li [p] ;
            y [0] -= Lx [p] * X [i] ;
            y [1] -= Lx [q] * X [i] ;
         }
         if(L->is_ll == true){
            y [0] /= d [0] ;
            y [1] = (y [1] - t * y [0]) / d [1] ;
         }
         else
            y [1] -= t * y [0] ;
         X [j  ] = y [0] ;
         X [j-1] = y [1] ;
         j -= 2 ;	    /* advance to the next column of L */

      }
      else
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of three columns of L */
         /* -------------------------------------------------------------- */

         double y [3], t [3] ;
         q = Lp [j-1] ;
         r = Lp [j-2] ;
         double d [3] ;
         d [0] = Lx [p] ;
         d [1] = Lx [q] ;
         d [2] = Lx [r] ;
         t [0] = Lx [q+1] ;
         t [1] = Lx [r+1] ;
         t [2] = Lx [r+2] ;
         if(L->is_ll == false){
            y [0] = X [j]   / d [0] ;
            y [1] = X [j-1] / d [1] ;
            y [2] = X [j-2] / d [2] ;
         }
         else{
            y [0] = X [j] ;
            y [1] = X [j-1] ;
            y [2] = X [j-2] ;
         }
         for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
         {
            int i = Li [p] ;
            y [0] -= Lx [p] * X [i] ;
            y [1] -= Lx [q] * X [i] ;
            y [2] -= Lx [r] * X [i] ;
         }
         if(L->is_ll == true){
            y [0] /= d [0] ;
            y [1] = (y [1] - t [0] * y [0]) / d [1] ;
            y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
         }
         else{
            y [1] -= t [0] * y [0] ;
            y [2] -= t [2] * y [0] + t [1] * y [1] ;
         }
         X [j-2] = y [2] ;
         X [j-1] = y [1] ;
         X [j  ] = y [0] ;
         j -= 3 ;	    /* advance to the next column of L */
      }
      //#endif
   }
   //t2 = clock();
   //clog<<"FBS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
}

// solve with ompenmp
// vector<Node_G *> children can be built when build_tree, so that
// memory will be more, bu the time for solving will be saved
//#if 0
void Circuit::solve_eq_pr(cholmod_factor *L, double *X){
   double *Lx;
   int *Li, *Lp, *Lnz;
   int p, q, r, lnz, pend;
   Lp = static_cast<int *>(L->p);
   Lx = static_cast<double*> (L->x);
   Li = static_cast<int*>(L->i) ;
   Lnz = static_cast<int *>(L->nz);
   int j, n = L->n ;
   // assign xp[i] = bnewp[i]
   //for(int i=0;i<n;i++)
	//X[i] = res[i]; 
   start_ptr_assign_1(); 

   clock_t t1, t2;
   t1 = clock();
   int count = 0;
   // FFS solve
   for (j = 0 ; j < n ; ){
      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      //if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
      {

         /* -------------------------------------------------------------- */
         /* solve with a single column of L */
         /* -------------------------------------------------------------- */
	      double y = X [j] ;
	      if(L->is_ll == true){
		      X[j] /= Lx [p] ;
	      }

	      if(X[j] != 0){
		 	count ++;
		      for (p++ ; p < pend ; p++)
		      {
			      X [Li [p]] -= Lx [p] * y ;
		      }

	      }
	      j++ ;	/* advance to next column of L */
      }
#if 0
      else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
      {

	      /* -------------------------------------------------------------- */
	      /* solve with a supernode of two columns of L */
	      /* -------------------------------------------------------------- */
	      //if(flag_col[j] == false && flag_col[j+1] == false){
		      double y [2] ;
		      q = Lp [j+1] ;
		      if(L->is_ll == true){
			      y [0] = X [j] / Lx [p] ;
			      y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
			      X [j  ] = y [0] ;
			      X [j+1] = y [1] ;
		      }

		      else{
			      y [0] = X [j] ;
			      y [1] = X [j+1] - Lx [p+1] * y [0] ;
			      X [j+1] = y [1] ;
		      }
		      for (p += 2, q++ ; p < pend ; p++, q++)
		      {
			      X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
		      }
	      //}
	      // else only  solve 1 column
	      //else if(flag_col[j] == false){

	      //}
	      j += 2 ;	    /* advance to next column of L */

      }
      else
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of three columns of L */
         /* -------------------------------------------------------------- */

         double y [3] ;
         q = Lp [j+1] ;
         r = Lp [j+2] ;
         if(L->is_ll == true){
            y [0] = X [j] / Lx [p] ;
            y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
            y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
            X [j  ] = y [0] ;
            X [j+1] = y [1] ;
            X [j+2] = y [2] ;
         }

         else{
            y [0] = X [j] ;
            y [1] = X [j+1] - Lx [p+1] * y [0] ;
            y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
            X [j+1] = y [1] ;
            X [j+2] = y [2] ;
         }
         for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
         {
            X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
         }
         j += 3 ;	    /* advance to next column of L */
      }
#endif
   }
   t2 = clock();
   clog<<count<<" out of: "<<n<<" are computed in FFS. "<<endl;
   clog<<"FFS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

   t1 = clock();
   // FBS solve
   for(j = n-1; j >= 0; ){

      /* get the start, end, and length of column j */
      p = Lp [j] ;
      lnz = Lnz [j] ;
      pend = p + lnz ;

      /* find a chain of supernodes (up to j, j-1, and j-2) */

      if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
      {

         /* -------------------------------------------------------------- */
         /* solve with a single column of L */
         /* -------------------------------------------------------------- */
         double d = Lx [p] ;
         if(L->is_ll == false)
            X[j] /= d ;
         for (p++ ; p < pend ; p++)
         {
            X[j] -= Lx [p] * X [Li [p]] ;
         }
         if(L->is_ll == true)
            X [j] /=  d ;
         j--;
      }
//#if 0
      else if (lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of two columns of L */
         /* -------------------------------------------------------------- */

         double y [2], t ;
         q = Lp [j-1] ;
         double d [2] ;
         d [0] = Lx [p] ;
         d [1] = Lx [q] ;
         t = Lx [q+1] ;
         if(L->is_ll == false){
            y [0] = X [j  ] / d [0] ;
            y [1] = X [j-1] / d [1] ;
         }
         else{
            y [0] = X [j  ] ;
            y [1] = X [j-1] ;
         }
         for (p++, q += 2 ; p < pend ; p++, q++)
         {
            int i = Li [p] ;
            y [0] -= Lx [p] * X [i] ;
            y [1] -= Lx [q] * X [i] ;
         }
         if(L->is_ll == true){
            y [0] /= d [0] ;
            y [1] = (y [1] - t * y [0]) / d [1] ;
         }
         else
            y [1] -= t * y [0] ;
         X [j  ] = y [0] ;
         X [j-1] = y [1] ;
         j -= 2 ;	    /* advance to the next column of L */

      }
      else
      {

         /* -------------------------------------------------------------- */
         /* solve with a supernode of three columns of L */
         /* -------------------------------------------------------------- */

         double y [3], t [3] ;
         q = Lp [j-1] ;
         r = Lp [j-2] ;
         double d [3] ;
         d [0] = Lx [p] ;
         d [1] = Lx [q] ;
         d [2] = Lx [r] ;
         t [0] = Lx [q+1] ;
         t [1] = Lx [r+1] ;
         t [2] = Lx [r+2] ;
         if(L->is_ll == false){
            y [0] = X [j]   / d [0] ;
            y [1] = X [j-1] / d [1] ;
            y [2] = X [j-2] / d [2] ;
         }
         else{
            y [0] = X [j] ;
            y [1] = X [j-1] ;
            y [2] = X [j-2] ;
         }
         for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
         {
            int i = Li [p] ;
            y [0] -= Lx [p] * X [i] ;
            y [1] -= Lx [q] * X [i] ;
            y [2] -= Lx [r] * X [i] ;
         }
         if(L->is_ll == true){
            y [0] /= d [0] ;
            y [1] = (y [1] - t [0] * y [0]) / d [1] ;
            y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
         }
         else{
            y [1] -= t [0] * y [0] ;
            y [2] -= t [2] * y [0] + t [1] * y [1] ;
         }
         X [j-2] = y [2] ;
         X [j-1] = y [1] ;
         X [j  ] = y [0] ;
         j -= 3 ;	    /* advance to the next column of L */
      }
//#endif
   }
   t2 = clock();
   //clog<<"FBS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
 }


// FFS solve
void Circuit::solve_col_FFS(cholmod_factor *L, double *X, int &j,
	double *Lx, int *Li, int *Lp, int *Lnz){
   int p, lnz, pend;
   //for(int i=0;i<L->n;i++) 
	//cout<<"i, X: "<<i<<" "<<X[i]<<endl;

   p = Lp [j] ;
   lnz = Lnz[j] ;
   pend = p + lnz ;

   double y = X [j] ;
   if(L->is_ll == true)
	   X[j] /= Lx [p] ;
   //cout<<"j, X: "<<j<<" "<<X[j]<<endl;

   for (++p ; p < pend ; p++){
	   X [Li [p]] -= Lx [p] * y ;
  }
}

// FBS solve
void Circuit::solve_col_FBS(cholmod_factor *L, double *X, int &j,
	double *Lx, int *Li, int *Lp, int *Lnz){
   int p, lnz, pend;
   
   p = Lp [j] ;
   lnz = Lnz [j] ;
   pend = p + lnz ;

   double d = Lx [p] ;
   if(L->is_ll == false)
	   X[j] /= d ;
   for (++p ; p < pend ; p++)
	   X[j] -= Lx [p] * X [Li [p]] ;
   if(L->is_ll == true)
	   X [j] /=  d ;
}

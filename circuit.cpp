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
#include "umfpack.h"
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"
using namespace std;

double Circuit::EPSILON = 1e-5;
double Circuit::OMEGA = 1.2;
double Circuit::OVERLAP_RATIO = 0.2;//0.2;
size_t Circuit::MAX_BLOCK_NODES = 4000;//100000;
int    Circuit::MODE = 0;
const int MAX_ITERATION = 1000;
const int SAMPLE_INTERVAL = 5;
const size_t SAMPLE_NUM_NODE = 10;

//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

vector<LAYER_DIR> Circuit::layer_dir(MAX_LAYER);

// constructor of Circuit class, name is optional
Circuit::Circuit(string _name):name(_name),
	x_min(INFTY),y_min(INFTY),x_max(0),y_max(0),
	circuit_type(UNKNOWN), VDD(0.0){
	// add ground node
	Node * gnd = new Node(string("0"), Point(-1,-1,-1));
	gnd->rep = gnd;
	this->add_node(gnd);

	for(int i=0;i<MAX_LAYER;i++)
		layer_dir[i]=NA;
	id_map = NULL;
	etree.clear();
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
	tree.clear();
	for(size_t i=0;i<etree.size();i++) 
		delete etree[i];
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
void Circuit::solve_init(Tran &tran){
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
	clog<<"number of nodes can be merged is: "<<count<<endl;
}
// partition the circuit to X_BLOCKS * Y_BLOCKS blocks
// according to the node size. Assuming they are distributed
// uniformly at random
void Circuit::partition_circuit(){
	size_t num_nodes = replist.size();
	size_t num_blocks =  num_nodes / MAX_BLOCK_NODES;
	if( num_nodes % MAX_BLOCK_NODES > 0 ) ++num_blocks;
	size_t len_x = x_max-x_min;
	size_t len_y = y_max-y_min;
	size_t X_BLOCKS, Y_BLOCKS;

	// Extreme case: only one node in x/y axis
	if( num_blocks == 1 ){
		X_BLOCKS = Y_BLOCKS = 1;
	}
	else if( len_x == 0 ){
		X_BLOCKS = 1;
		Y_BLOCKS = num_blocks;
	}
	else if (len_y == 0 ){
		Y_BLOCKS = 1;
		X_BLOCKS = num_blocks;
	}
	else{// compute length of x and y according to their ratio
		double ratio = double(len_x) / double(len_y);
		X_BLOCKS = sqrt(num_blocks/ratio);
		Y_BLOCKS = X_BLOCKS*ratio;
		if( X_BLOCKS * Y_BLOCKS < num_blocks ) Y_BLOCKS+=1; 
		num_blocks = X_BLOCKS * Y_BLOCKS;
	}
	//X_BLOCKS =2;
	//Y_BLOCKS =2;
	clog<<"num_nodes: "<<num_nodes<<endl;
	clog<<"num_blocks: "<<X_BLOCKS<<" / "<<Y_BLOCKS <<endl;
	block_info.X_BLOCKS = X_BLOCKS;
	block_info.Y_BLOCKS = Y_BLOCKS;
	block_info.resize(X_BLOCKS * Y_BLOCKS);
}

// build up block info
// 1. Find block divide point
// 2. Set each node in replist into block
// 3. Compute block size
// 4. Insert boundary netlist into map
void Circuit::block_init(){
	block_info.set_len_per_block(x_min, x_max, y_min, y_max, OVERLAP_RATIO);
	block_info.update_block_geometry();
	find_block_size();
	copy_node_voltages_block();
}

void Circuit::block_boundary_insert_net(Net * net){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	if(nd[0]->is_ground() || nd[1]->is_ground()) return;

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};

	// idea: for each block where node_a is in
	// find whether node_b is also in the block
	// if yes, then this net is not a boundary net
	// for this block
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		// for each block k in ls[j]
		for(size_t k=0;k<p->size();k++){
			// find whether block_id in another list
			size_t block_id = (*p)[k];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				Block & blk = block_info[block_id];
				blk.boundary_netlist.push_back(net);
			}
		}// end of for k
	}// end of for j
}

// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix(Matrix *A){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			for(size_t i=0;i<ns.size();i++){
				assert( fzero(ns[i]->value) == false );
				stamp_block_resistor(ns[i], A);
			}
			break;
		case CURRENT:
			for(size_t i=0;i<ns.size();i++){
				stamp_block_current(ns[i]);
			}
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD(ns[i], A);
			}
			break;
		case CAPACITANCE:	
			break;
		case INDUCTANCE:
			// rephase inductance net
			for(size_t i=0;i<ns.size();i++){
				stamp_block_inductance_dc(ns[i], A);
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
        make_A_symmetric_block();

	// after stamping, convert A to column compressed form
	for(size_t i=0;i<num_blocks;i++){
		if(block_info[i].count>0){
		        //clog<<endl<<A[i]<<endl;
			A[i].set_row(block_info[i].count);
			block_info[i].CK_decomp(A[i], cm);
		}
	}
}

// 1. mark rep nodes into corresponding blocks
// 2. find block size of replist nodes
// 3. allocate rhs size of each block
// 4. find local block index for each node
void Circuit::find_block_size(){
	// start assign replist node into blocks
	size_t n = replist.size();
	Node * p = NULL;
	
	// find which blocks the node belongs to,
	// and its index in the block
	for(size_t i=0;i<n;i++){
		p=replist[i];
		set_blocklist(p);	
	}

	// for each block, allocate resource
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		block.allocate_resource(cm);
	}
}

void Circuit::solve(Tran &tran){
	solve_LU(tran);
	clog<<"finish solve. "<<endl;
}

// solve Circuit
bool Circuit::solve_IT(Tran &tran){
	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	solve_init(tran);
	/*if( replist.size() <= 2*MAX_BLOCK_NODES ){
		clog<<"Use direct LU instead."<<endl;
		solve_LU_core();
		return true;
	}*/

	//select_omega();
	partition_circuit();

	// Only one block, use direct LU instead
	/*if( block_info.size() == 1 ){
		clog<<"Use direct LU instead."<<endl;
		solve_LU_core(tran);
		return true;
	}*/

        cm = &c;
        cholmod_start(cm);
        cm->print = 5;

	block_init();
	
	clog<<"e="<<EPSILON
	    <<"\to="<<OMEGA
	    <<"\tr="<<OVERLAP_RATIO
	    <<"\tb="<<MAX_BLOCK_NODES
	    <<"\tmode="<<MODE<<endl;

	// first solve DC part
	num_blocks = block_info.size();
	Matrix A[num_blocks];
      
	bool successful = false;	
	successful = solve_IT_dc(A);
	//cout<<"dc solution for ckt:  "<<get_name()<<endl;
	//cout<<nodelist<<endl;	
	//return successful;
	
	// clear all elements in b for tr restamp
	for(size_t i=0;i<num_blocks;i++){
	    for(size_t j=0;j<block_info[i].count;j++)	
	       block_info[i].bp[j]=0;
	}
	successful = solve_IT_tr(A, tran);

	return successful;
}

// TODO: add comment
void Circuit::node_voltage_init(){
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		for(size_t j=0;j<block.count;j++){
			block.xp[i] = VDD;
			block.nodes[i]->value = VDD;
		}
	}
}

// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration(){	
   double diff = .0, max_diff = .0;
   int flag_tr = 0; // for dc, flag_tr = 0;
   for(size_t i=0;i<block_info.size();i++){
      Block &block = block_info[i];
      if( block.count == 0 ) continue;

      // b_new is the final rhs
      block.update_rhs(flag_tr);

      // backup the old voltage value
      block.x_old = block.xp;
      block.solve_CK(cm);
      block.xp = static_cast<double*>(block.x_ck->x);

      // modify node voltage with OMEGA and old voltage value
      diff = modify_voltage(block, block.x_old);

      if( max_diff < diff ) max_diff = diff;
   }
   return max_diff;
}

double Circuit::solve_iteration_tr(Tran &tran, double &time){
	// stamp into b_new
	stamp_block_current_tr(time);

	// stamp into b_new
	modify_block_rhs_tr(tran);

	int counter = 0;
	double diff=0;
	bool successful = false;

	while( counter++ < MAX_ITERATION ){
		diff = solve_iteration_tr_step();

		//clog<<"iter, diff: "<<counter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}

	get_voltages_from_block_LU_sol();
	return successful;
}

double Circuit::modify_voltage(Block & block, double * x_old){
	double max_diff = 0.0;
	OMEGA = 1.0;
	for(size_t i=0;i<block.count;i++){
		block.xp[i] = (1-OMEGA) * x_old[i] + OMEGA * block.xp[i];
		block.nodes[i]->value = block.xp[i];
		double diff = fabs(x_old[i] - block.xp[i]);
		if( diff > max_diff ) max_diff = diff;
	}
	return max_diff;
}

// TODO: add comment
double Circuit::modify_voltage(Block & block){
	double max_diff = 0.0;
	for(size_t i=0;i<block.count;i++){
		Node * nd = block.nodes[i];
		double old = nd->value;

		block.xp[i] = (1-OMEGA) * old + OMEGA * block.xp[i];

		// * copy back to the actual node *
		nd->value = block.xp[i];
		double diff = fabs(old-block.xp[i]);
		if( diff > max_diff ) max_diff = diff;
	}
	return max_diff;
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
   //return;

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
   
   stamp_current_tr(bp, tran, time);
   
   double t1, t2;
   t1 = omp_get_wtime();
   Algebra::CK_decomp(A, L, cm);
   Lp = static_cast<int *>(L->p);
   Lx = static_cast<double*> (L->x);
   Li = static_cast<int*>(L->i) ;
   Lnz = static_cast<int *>(L->nz);

   //cholmod_print_factor(L, "L", cm);
   
   t2 = omp_get_wtime();
   clog<<"decomp cost: "<<1.0*(t2-t1)<<endl;
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
   modify_rhs_tr(bnewp, xp, tran, iter);

   t1 = omp_get_wtime(); 
   solve_eq(xp); 
   t2 = omp_get_wtime();
   clog<<"time for solve_tr: "<<t2-t1<<endl;

   save_tr_nodes(tran, xp);
   time += tran.step_t;
   t1 = omp_get_wtime();
   // then start other iterations
   while(time < tran.tot_t){// && iter < 2){
       // bnewp[i] = bp[i];
       for(int i=0;i<n;i++)
	bnewp[i] = bp[i];
	//start_ptr_assign();


      // only stamps if net current changes
      // set bp into last state
      //stamp_current_tr(bnewp, tran, time);
      stamp_current_tr_1(bp, bnewp, tran, time);
     // get the new bnewp
      modify_rhs_tr(bnewp, xp, tran, iter);
    
	
      solve_eq(xp);
      //x = cholmod_solve(CHOLMOD_A, L, bnew, cm); 

      save_tr_nodes(tran, xp);
      time += tran.step_t;
      iter ++;
   }
   t2 = omp_get_wtime();
   clog<<"1000 iters cost: "<<t2-t1<<endl;

   release_tr_nodes(tran);
   cholmod_free_dense(&x, cm);
   cholmod_free_dense(&b, cm);
   cholmod_finish(&c);
}

// solve the node voltages using direct LU
void Circuit::solve_LU(Tran &tran){
	solve_init(tran);
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

// copy solution of block into circuit
void Circuit::get_voltages_from_block_LU_sol(){

	size_t block_id;
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		block_id = node->rep->blocklist[0];
		Block &block = block_info[block_id];
		size_t id = node->rep->id_in_block[0];
		double v = block.xp[id];		// get its rep's value
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

// 1. copy node voltages from the circuit to a Vec
//    from = true then copy circuit to x
//    else copy from x to circuit
// 2. map block voltage into global_index
void Circuit::copy_node_voltages_block(bool from){
	size_t id;
	double value = 0;
	if( from == true ){
		for(size_t i=0;i<replist.size();i++){
			Node *node = replist[i];
			const vector<size_t> &block_id = node->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->id_in_block[j];
				value = replist[i]->value;
				if(replist[i]->nbr[TOP]!=NULL 
				  && replist[i]->nbr[TOP]->type
				  == INDUCTANCE){
				   value = replist[i]->nbr[TOP]->ab[0]->value;
				   if(replist[i]->nbr[TOP]->ab[0]->isS()!=Y)
				      value = replist[i]->nbr[TOP]->ab[1]->value;
                                }

				block.xp[id] = value;
				//clog<<"rep, id, xp: "<<*replist[i]<<" "<<id<<" "<<block.xp[id]<<endl;
				block.nodes[id] = replist[i];
			}
		}
	}
	else{
		for(size_t i=0;i<nodelist.size()-1;i++){
			Node * node = nodelist[i];
			const vector<size_t> &block_id = 
				node->rep->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->rep->id_in_block[j];
				node->value = block.xp[id];
			}
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
void Circuit::stamp_current_tr(double *b, Tran &tran, double &time){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type == CURRENT)
			for(size_t i=0;i<ns.size();i++)
				stamp_current_tr_net(b, ns[i], time);
	}
}

// stamp transient current values into rhs
void Circuit::stamp_current_tr_1(double *bp, double *b, Tran &tran, double &time){
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
void Circuit::modify_rhs_tr(double * b, double *x, Tran &tran, int &iter){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++)
				modify_rhs_c_tr(ns[i], b, x, tran, iter);
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_rhs_l_tr(ns[i], b, x, tran, iter);	
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
void Circuit::modify_rhs_c_tr(Net *net, double * rhs, double *x, Tran &tran, int &iter){
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
void Circuit::modify_rhs_l_tr(Net *net, double *rhs, double *x, Tran &tran, int &iter){
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

// =========== stamp block version of matrix =======

void Circuit::stamp_block_resistor(Net * net, Matrix * A){
   int boundary_flag = 0;
   Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
   if(nd[0]->is_ground() || nd[1]->is_ground()) 
      boundary_flag = 1;;

   double G;	
   G = 1./net->value;
   //clog<<"net: "<<*net<<endl;

   vector<size_t> * ls[] = {&nd[0]->blocklist, 
      &nd[1]->blocklist};
   vector<size_t>::const_iterator it;
   for(size_t j=0;j<2;j++){
      vector<size_t> *p = ls[j], *q = ls[1-j];
      Node *nk = nd[j], *nl = nd[1-j];
      //clog<<*nk<<" "<<*nl<<endl;
      for(size_t i=0;i<p->size();i++){
         // find whether block_id in another list
         size_t block_id = (*p)[i];
         it = find( (*q).begin(), (*q).end(), block_id);

         // block_id not in another list, net is boundary
         if( it == (*q).end() ){
            Block & blk = block_info[block_id];
            if(boundary_flag ==0)
               blk.boundary_netlist.push_back(net);
            if( !nk->is_ground() && nk->isS() != Y 
              &&(nk->nbr[TOP]==NULL || 
                 nk->nbr[TOP]->type !=INDUCTANCE) ) {
               // stamp value into block_ids
               size_t k1 = nk->id_in_block[i];
               Matrix &pk = A[block_id];	
               pk.push_back(k1,k1, G);
               //clog<<k1<<" "<<k1<<" "<<G<<endl;
            }
         }
         // else 2 nodes belongs to the same block
         // stamp resistor
         else if( !nk->is_ground() && nk->isS()!=Y &&
           (nk->nbr[TOP]==NULL || 
            nk->nbr[TOP]->type !=
            INDUCTANCE)) {
            size_t k1 = nk->id_in_block[i];
            Matrix &pk = A[block_id];

            size_t j1 = it - (*q).begin();
            size_t l1 = nl->id_in_block[j1];

            pk.push_back(k1,k1, G);
            //clog<<k1<<" "<<k1<<" "<<G<<endl;
            if(!nl->is_ground() &&(nl->nbr[TOP]==NULL || 
                 nl->nbr[TOP]->type != INDUCTANCE)&&(l1 > k1)){
               pk.push_back(k1,l1,-G);
               //clog<<k1<<" "<<l1<<" "<<-G<<endl;
            }
         }
      }// end of for k
   }// end of for j	
}

// node in resistor net can be gound node here
void Circuit::stamp_block_resistor_tr(Net * net, Matrix * A){
   Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
   //if(nd[0]->is_ground() || nd[1]->is_ground()) return;

   double G;	
   G = 1./net->value;

   vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
   vector<size_t>::const_iterator it;
   for(size_t j=0;j<2;j++){
      vector<size_t> *p = ls[j], *q = ls[1-j];
      Node *nk = nd[j], *nl = nd[1-j];
      for(size_t i=0;i<p->size();i++){
         // find whether block_id in another list
         size_t block_id = (*p)[i];
         it = find( (*q).begin(), (*q).end(), block_id);

         // block_id not in another list, net is boundary
         if( it == (*q).end() ){
            if( nk->isS() != Y && !nk->is_ground()//){
               &&nk->nbr[TOP]!=NULL &&
                 nk->nbr[TOP]->type ==INDUCTANCE) {
                    // stamp value into block_ids
                    size_t k1 = nk->id_in_block[i];
                    Matrix &pk = A[block_id];	
                    pk.push_back(k1,k1, G);
                 }
         }
         // else 2 nodes belongs to the same block
         // stamp resistor
            else if( nk->isS()!=Y && !nk->is_ground()//){
               &&nk->nbr[TOP]!=NULL &&
                 nk->nbr[TOP]->type ==INDUCTANCE) {
                    size_t k1 = nk->id_in_block[i];
                    Matrix &pk = A[block_id];

                    size_t j1 = it - (*q).begin();
                    size_t l1 = nl->id_in_block[j1];

                    pk.push_back(k1,k1, G);
                    if(!nl->is_ground()&&
                      (nl->nbr[TOP]==NULL || 
                       nl->nbr[TOP]->type != INDUCTANCE)){
                       if(l1 > k1)
                          pk.push_back(k1,l1,-G);
                       else if(l1 < k1)
                          pk.push_back(l1,k1,-G);
                    }
                 }
      }// end of for k
   }// end of for j	
}

// stamp inductance Geq = delta_t/(2L)
// stamp capacitance Geq = 2C / delta_t
// node in resistor net can be gound node here
void Circuit::stamp_block_resistor_cl_tr(Matrix *A, Net * net, double &Geq){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		Node *nk = nd[j], *nl = nd[1-j];
		for(size_t i=0;i<p->size();i++){
			// find whether block_id in another list
			size_t block_id = (*p)[i];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				//blk.boundary_netlist.push_back(net);
				if( nk->isS() != Y && !nk->is_ground()) {
					// stamp value into block_ids
					size_t k1 = nk->id_in_block[i];
					Matrix &pk = A[block_id];	
					pk.push_back(k1,k1, Geq);
					if(net->type == INDUCTANCE)
						pk.push_back(k1,k1,-1);
				}
			}
			// else 2 nodes belongs to the same block
			// stamp resistor
			else if( nk->isS()!=Y && !nk->is_ground()) {
				size_t k1 = nk->id_in_block[i];
				Matrix &pk = A[block_id];

				size_t j1 = it - (*q).begin();
				size_t l1 = nl->id_in_block[j1];

				pk.push_back(k1,k1, Geq);
				if(net->type == INDUCTANCE)
					pk.push_back(k1,k1,-1);

				if(nl->isS()!=Y && !nl->is_ground() && k1 <l1)
					pk.push_back(k1,l1,-Geq);
			}
		}// end of for k
	}// end of for j	
}

void Circuit::stamp_block_capacitance_dc(Net * net, Matrix * A){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j];
		Node *nk = nd[j], *nl = nd[1-j];
		for(size_t i=0;i<p->size();i++){
			// find whether block_id in another list
			size_t block_id = (*p)[i];
			// 2 nodes belongs to the same block
			// stamp short net
			if( nk->isS()!=Y && !nk->is_ground()) {
				size_t k1 = nk->id_in_block[i];
				Matrix &pk = A[block_id];
				size_t l1 = nl->id_in_block[i];

				pk.push_back(k1,k1, 0);
				if(!nl->is_ground())
					pk.push_back(k1, l1, 0);
			}
		}// end of for k
	}// end of for j	
}

void Circuit::stamp_block_inductance_dc(Net * net, Matrix * A){
   //clog<<"net: "<<*net<<endl;	
      Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j];
		Node *nk = nd[j], *nl = nd[1-j];
		for(size_t i=0;i<p->size();i++){
			// find whether block_id in another list
			size_t block_id = (*p)[i];
			// 2 nodes belongs to the same block
			// stamp short net
			if( nk->isS()!=Y && !nk->is_ground()) {
				size_t k1 = nk->id_in_block[i];
				Matrix &pk = A[block_id];
				size_t l1 = nl->id_in_block[i];

				pk.push_back(k1,k1, 1);
				if(!nl->is_ground()){
				// general stamping
				// make it symmetric
				 block_info[block_id].bp[k1] = 
				   block_info[block_id].bp[l1];
				 //clog<<block_info[block_id].bp[k1]<<" "<<block_info[block_id].bp[l1]<<endl;
                                }
			}
		}// end of for k
	}// end of for j	
}

void Circuit::stamp_block_current(Net * net){
   //clog<<"net: "<<*net<<endl;	
   Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && nk->isS()!=Y ){//&&
		//(nk->nbr[TOP]==NULL || 
	 	//(nk->nbr[TOP]->type !=INDUCTANCE))) { 
		for(size_t i=0;i<nk->blocklist.size();i++){
			size_t block_idk = nk->blocklist[i];
			Block &block_k = block_info[block_idk];
			size_t k = nk->id_in_block[i];
			block_k.bp[k] += -net->value;
			//clog<<"insert: "<<k<<" "<<block_k.bp[k]<<endl;
		}
	}
	if( !nl->is_ground() && nl->isS()!=Y ){//&&
		//(nl->nbr[TOP]==NULL || 
		//(nl->nbr[TOP]->type !=INDUCTANCE))) {
		for(size_t i=0;i<nl->blocklist.size();i++){
			size_t block_idl = nl->blocklist[i];
			Block & block_l = block_info[block_idl];
			size_t l = nl->id_in_block[i];
			block_l.bp[l] +=  net->value;
			//clog<<"insert: "<<l<<" "<<block_l.bp[l]<<endl;
		}
	}
}

// stamp Ieq for capcitance and inductance
void Circuit::stamp_block_current_tr(double &time){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type == CURRENT)
			for(size_t i=0;i<ns.size();i++)
				stamp_block_current_tr_net(ns[i], time);
	}
}

void Circuit::stamp_block_current_tr_net(Net *net, double &time){
	current_tr(net, time);
	//clog<<"net: "<<*net<<endl;
	//clog<<"current: "<<current<<endl;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && nk->isS()!=Y) { 
		for(size_t i=0;i<nk->blocklist.size();i++){
			size_t block_idk = nk->blocklist[i];
			//clog<<"block_idk: "<<block_idk<<endl;
			Block &block_k = block_info[block_idk];
			size_t k = nk->id_in_block[i];
			block_k.b_trp[k] += -net->value;//current;
			//clog<<"k, b_tr: "<<k<<" "<<block_k.b_trp[k]<<endl;
		}
	}
	if( !nl->is_ground() && nl->isS()!=Y) {
		for(size_t i=0;i<nl->blocklist.size();i++){
			size_t block_idl = nl->blocklist[i];
			//clog<<"block_idl: "<<block_idl<<endl;
			Block & block_l = block_info[block_idl];
			size_t l = nl->id_in_block[i];
			block_l.b_trp[l] +=  net->value;//current;
			//clog<<"l, b_tr: "<<l<<" "<<block_l.b_trp[l]<<endl;
		}
	}
}

// add Ieq into rhs
// Ieq = i(t)+2C/delta_t *v(t)
void Circuit::modify_block_rhs_c_tr(Net *net, Tran &tran){
	double i_t = 0;
	double temp = 0;
	double Ieq = 0;
	
	//clog<<"c net: "<<*net<<endl;
	
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	if(nk->isS() != Z)
		swap<Node *>(nk, nl);

	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;

	Net *r = nk->nbr[TOP];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to Z node
	if(a->isS()!=Z) swap<Node *>(a, b);
	//clog<<"a, b: "<<*a<<" "<<*b<<endl;

	i_t = (b->value - a->value) / r->value;
	//clog<<"i_t: "<<i_t<<endl;
	temp = 2*net->value / tran.step_t * 
		(nk->value - nl->value);
	//clog<<"nk-nl "<<(nk->value - nl->value)<<" "<<2*net->value/tran.step_t<<" "<<temp<<endl;
	
	Ieq  = (i_t + temp);
	//clog<< "Ieq is: "<<Ieq<<endl;
	if(!nk->is_ground()&& nk->isS()!=Y){
		for(size_t i=0;i<nk->blocklist.size();i++){
			size_t block_idk = nk->blocklist[i];
			//clog<<"block_idk: "<<block_idk<<endl;
			Block &block_k = block_info[block_idk];
			size_t k = nk->id_in_block[i];
			block_k.b_trp[k] += Ieq;
			//clog<<*nk<<" rhs +: "<<block_k.b_trp[k]<<endl;
		}

	}
	if(!nl->is_ground()&& nl->isS()!=Y){
		for(size_t i=0;i<nl->blocklist.size();i++){
			size_t block_idl = nl->blocklist[i];
			//clog<<"block_idl: "<<block_idl<<endl;
			Block & block_l = block_info[block_idl];
			size_t l = nl->id_in_block[i];
			block_l.b_trp[l] +=  -Ieq;
			//clog<<*nl<<" rhs +: "<<block_l.b_trp[l]<<endl;
		}
	}
}

// add Ieq into rhs
// Ieq = i(t) + delta_t / (2*L) *v(t)
void Circuit::modify_block_rhs_l_tr(Net *net, Tran &tran){
	//clog<<"l net: "<<*net<<endl;
	Node *nk = net->ab[0]->rep;
	Node *nl = net->ab[1]->rep;
	// nk point to X node
	if(nk->isS() !=X) 
		swap<Node*>(nk, nl);
	//clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
	double Ieq = 0;

	double i_t = 0;
	double temp = 0;
	temp = tran.step_t / (2*net->value) * 
		(nl->value - nk->value);
	//clog<<"delta_t/2L, nl-nk, temp: "<<tran.step_t / (2*net->value)<<" "<<(nl->value-nk->value)<<" "<<temp<<endl;
	
	Net *r = nk->nbr[BOTTOM];
	Node *a = r->ab[0]->rep;
	Node *b = r->ab[1]->rep;
	// a point to X node
	if(a->isS()!=X) swap<Node*>(a, b);
	i_t = (a->value - b->value) / r->value;
	//clog<<"resiste r: "<<*r<<endl;
	//clog<<*a<<" "<<*b<<endl;
	//clog<<"a, b, r, i_t: "<<a->value<<" "<<b->value<<" "<<
		//r->value<<" "<<i_t<<endl;
	Ieq  = i_t + temp;
	//clog<<"Ieq: "<<Ieq<<endl;
	if(nk->isS() !=Y && !nk->is_ground()){
		for(size_t i=0;i<nk->blocklist.size();i++){
			size_t block_idk = nk->blocklist[i];
			//clog<<"block_idk: "<<block_idk<<endl;
			Block &block_k = block_info[block_idk];
			size_t k = nk->id_in_block[i];
			//clog<<"trp bef: "<<block_k.b_trp[k]<<endl;
			block_k.b_trp[k] += Ieq;
			//clog<<*nk<<" rhs +: "<<block_k.b_trp[k]<<endl;
		}
	}
	if(nl->isS()!=Y && !nl->is_ground()){
		for(size_t i=0;i<nl->blocklist.size();i++){
			size_t block_idl = nl->blocklist[i];
			//clog<<"block_idl: "<<block_idl<<endl;
			Block & block_l = block_info[block_idl];
			size_t l = nl->id_in_block[i];
			block_l.b_trp[l] +=  -Ieq;
			//clog<<*nl<<" rhs-: "<<block_l.b_trp[l]<<endl;
		}
	}
}

void Circuit::stamp_block_VDD(Net * net, Matrix * A){
	// find the non-ground node
	Node * X = net->ab[0];
	//clog<<"net: "<<*net<<endl;

	if( X->is_ground() ) X = net->ab[1];
	for(size_t i=0;i<X->rep->id_in_block.size();i++){
		size_t block_id = X->rep->blocklist[i];
		Block &block = block_info[block_id];	
		Matrix & p = A[block_id];
		size_t id =X->rep->id_in_block[i];
		p.push_back(id, id, 1.0);
		//clog<<id<<" "<<id<<" "<<1.0<<endl;
		Net * south = X->rep->nbr[SOUTH];
		if( south != NULL &&
	   	 south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
			assert( feqn(1.0, block.bp[id]) ); 
			block.bp[id] = net->value;	    // modify it
		}
		else{
			block.bp[id] += net->value;
		}
	}
}

void Circuit::stamp_block_VDD_tr(Net * net){
	// find the non-ground node
	Node * X = net->ab[0];

	if( X->is_ground() ) X = net->ab[1];
	for(size_t i=0;i<X->rep->id_in_block.size();i++){
		size_t block_id = X->rep->blocklist[i];
		Block &block = block_info[block_id];	
		//Matrix & p = A[block_id];
		size_t id =X->rep->id_in_block[i];
		//p.push_back(id, id, 1.0);
		Net * south = X->rep->nbr[SOUTH];
		if( south != NULL &&
	   	 south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
			assert( feqn(1.0, block.bp[id]) ); 
			block.bp[id] = net->value;	    // modify it
		}
		else{
			block.bp[id] += net->value;
		}
	}
}


// set the block_id in a node
// according to the following order
// where 0 is the original block the node should be 
// 8 1 2
// 7 0 3
// 6 5 4
//
void Circuit::set_blocklist(Node * nd){
	const double len_per_block_x = block_info.len_per_block_x;
	const double len_per_block_y = block_info.len_per_block_y;
	const double len_ovr_x = block_info.len_ovr_x;
	const double len_ovr_y = block_info.len_ovr_y;
	const size_t X_BLOCKS = block_info.X_BLOCKS;
	const size_t Y_BLOCKS = block_info.Y_BLOCKS;
	const long x = nd->pt.x - x_min;	// point relative location
	const long y = nd->pt.y - y_min;
	size_t bx0 = x / len_per_block_x;
	size_t by0 = y / len_per_block_y;
	long bx, by;		// block index
	double lx, ly, ux, uy;	// block bounding box

	const long dx[]={0, 0, 1, 1,  1,  0, -1, -1, -1};
	const long dy[]={0, 1, 1, 0, -1, -1, -1,  0, 1};

	// test whether the node is in one of the nine blocks
	for(int i=0;i<9;i++){
		bx = bx0 + dx[i];
		by = by0 + dy[i];

		// check the block index is valid
		if( bx < 0 || bx >= (long)X_BLOCKS ) continue;
		if( by < 0 || by >= (long)Y_BLOCKS ) continue;

		// compute block coordinate
		lx = bx * len_per_block_x - len_ovr_x;
		ly = by * len_per_block_y - len_ovr_y;
		ux = (bx+1) * len_per_block_x + len_ovr_x;
		uy = (by+1) * len_per_block_y + len_ovr_y;

		// check if the point is in the block
		if( !(x>=lx && x<=ux && y>=ly && y<=uy) ) continue;	

		size_t id = by * X_BLOCKS + bx;
		assert( id<X_BLOCKS* Y_BLOCKS );

		/*
		vector<size_t >::const_iterator it;
		it = find(nd->blocklist.begin(),nd->blocklist.end(), id);
		if(it!=nd->blocklist.end()){
			printf("id=%ld, i=%d\n", id,i);
			printf("xbase,ybase=%lf %lf\n", lx_base, ly_base);
			printf("lx,ly=%lf %lf\n", lx, ly);
			printf("ux,uy=%lf %lf\n", ux, uy);
			printf("x,y=%ld %ld\n", x, y);
			printf("bx,by=%ld %ld\n", bx, by);
			printf("xc,yc=%ld %ld\n", x_center, y_center);
			continue;
		}
		*/
		nd->blocklist.push_back(id);
		Block &block = block_info[id];
		nd->id_in_block.push_back(block.count++);
	}
}

void Circuit::get_parameters(
		double & epsilon,
		double & omega,
		double & overlap_ratio,
		size_t & max_block_nodes,
		int & mode){
	epsilon		= EPSILON;
	omega		= OMEGA;
	overlap_ratio	= OVERLAP_RATIO; 
	max_block_nodes	= MAX_BLOCK_NODES;
	mode		= MODE;
}

// default values of these parameters are at the begining of this file
void Circuit::set_parameters(
		double epsilon, 
		double omega, 
		double overlap_ratio,
		size_t max_block_nodes,
		int mode){
	EPSILON		= epsilon;
	OMEGA		= omega;
	OVERLAP_RATIO 	= overlap_ratio;
	MAX_BLOCK_NODES	= max_block_nodes;
	MODE		= mode;
}

// choose an appropriate omega for the circuit s.t.
// - node size (use replist)
// - type (c4 or wb)
// - number of layers
void Circuit::select_omega(){
	double omega=OMEGA;
	size_t num_nodes = replist.size();
	size_t num_layers = layers.size();
	if( num_nodes < 0.05e6 )
		omega=1.0;
	else if (num_nodes < 0.2e6 )
		omega = 1.1;
	else if (num_nodes < 0.3e6 )
		omega = 1.2;
	else if (num_nodes < 0.5e6 )
		omega = 1.3;
	else if (num_nodes < 1.2e6 )
		omega = 1.4;
	else
		omega = 1.5;

	if( circuit_type == WB && num_layers >= 8 ) omega += 0.2;

	if( circuit_type == C4 ) omega += 0.1;

	if( name == "GND" && num_nodes < 1.2e6) omega -= 0.1;

	if( omega >= 1.6 ) omega = 1.6;
	if( omega <= 1.0 ) omega = 1.0;

	OMEGA = omega;
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

bool Circuit:: solve_IT_dc(Matrix *A){
   // first solve DC part
   stamp_block_matrix(A);

   int counter = 0;
   double diff=0;
   bool successful = false;

   while( counter++< MAX_ITERATION ){
      diff = solve_iteration();

      //clog<<"iter, diff: "<<counter<<" "<<diff<<endl;
      if( diff < EPSILON ){
         successful = true;
         break;
      }
   }

   clog<<"# iter: "<<counter<<endl;
   get_voltages_from_block_LU_sol();
   return successful;
}

bool Circuit:: solve_IT_tr(Matrix *A, Tran &tran){
	// stamp again for transient analysis
	stamp_block_matrix_tr(A, tran);

	double diff=0;
	bool successful = false;

	double time = 0;
	double *x=NULL;
	int iter = 0;
	while(time < tran.tot_t ){//&& iter<1){
		for(size_t i=0;i<num_blocks;i++){ 
                   for(size_t k=0;k<block_info[i].count;k++){
			block_info[i].b_trp[k] = block_info[i].bp[k];
		        //if(i==0)
		        //clog<<"bp first: "<<block_info[i].b_trp[k]<<endl;
                   }
		}

                make_A_symmetric_block_tr(tran);
                   //for(size_t k=0;k<block_info[0].count;k++)
                      //clog<<"b_trp after sym: "<<block_info[0].b_trp[k]<<endl;
		diff = solve_iteration_tr(tran, time);
		//cout<<"time: "<<time<<endl;
		save_tr_nodes(tran, x);
		time += tran.step_t;
		iter++;
	}
	//print_tr_nodes(tran);
	return successful;
}

// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix_tr(Matrix *A, Tran &tran){
	double Geq = 0;
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		switch(type){
		case RESISTOR:
			for(size_t i=0;i<ns.size();i++){
				assert( fzero(ns[i]->value) == false );
				stamp_block_resistor_tr(ns[i], A);
			}
			break;
		case CURRENT:
			//for(size_t i=0;i<ns.size();i++){
				//stamp_block_current(ns[i], A);
			//}
			break;
		case VOLTAGE:
			for(size_t i=0;i<ns.size();i++){
				if( fzero(ns[i]->value)  && 
				    !ns[i]->ab[0]->is_ground() &&
				    !ns[i]->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD_tr(ns[i]);
			}
			break;
		case CAPACITANCE: // Geq = 2C/delta_t
			for(size_t i=0;i<ns.size();i++){
				Geq = 2*ns[i]->value / tran.step_t;
				stamp_block_resistor_cl_tr(A, 
						ns[i], Geq);
			}
			break;
		case INDUCTANCE: // Geq = delta_t /(2L)
			for(size_t i=0;i<ns.size();i++){
				Geq = tran.step_t/(2*ns[i]->value);
				stamp_block_resistor_cl_tr(A, 
						ns[i], Geq);
			}
			break;

		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}

	// after stamping, convert A to column compressed form
	for(size_t i=0;i<num_blocks;i++){
		if(block_info[i].count>0){
			//A[i].set_row(block_info[i].count);
			block_info[i].CK_decomp(A[i], cm);
		}
	}
}

// update block rhs by transient nets
void Circuit::modify_block_rhs_tr(Tran &tran){
	//for(size_t i=0;i<num_blocks;i++)
		//clog<<"b_tr: "<<block_info[i].b_tr<<endl;
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetPtrVector & ns = net_set[type];
		if(type ==CAPACITANCE){	
			for(size_t i=0;i<ns.size();i++)
				modify_block_rhs_c_tr(ns[i], 
						tran);
			//for(int i=0;i<block_info[0].count;i++)
			   //clog<<"b_trp after cap: "<<block_info[0].b_trp[i]<<endl;
		}
		else if(type == INDUCTANCE){
			for(size_t i=0;i<ns.size();i++){
				modify_block_rhs_l_tr(ns[i],
						tran);	
			}
                        //for(int i=0;i<block_info[0].count;i++)
			   //clog<<"b_trp after ind: "<<block_info[0].b_trp[i]<<endl;

		}
	}
}

double Circuit::solve_iteration_tr_step(){
	double diff = .0, max_diff = .0;
	int flag_tr =1;
	//clog<<endl<<endl<<endl;
	for(size_t i=0;i<block_info.size();i++){
		Block &block = block_info[i];
		if( block.count == 0 ) continue;
		
		// b_tr is the stamped rhs of each time step
		block.update_rhs(flag_tr);

		// backup the old voltage value
		for(size_t k=0;k<block.count;k++)
		  block.x_old[k] = block.xp[k];
		block.solve_CK(cm);
		block.xp = static_cast<double*>(block.x_ck->x);

		// modify node voltage with OMEGA and old voltage value
		diff = modify_voltage(block, block.x_old);

		if( max_diff < diff ) max_diff = diff;
	}
	return max_diff;
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

void Circuit::make_A_symmetric_block(){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	vector<size_t> *p, *q;
	Node *nk, *nl, *nr;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
	
		Net *net = *it;
		Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};	
		//clog<<"symmetric net: "<<*net<<endl;
		if(nd[0]->isS()!=X && nd[1]->isS()!=X) continue;
		double G;
		G = 1./net->value;
		//clog<<"symmetric net again: "<<*net<<endl;

		vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
		vector<size_t >::const_iterator it_1;
		for(size_t j=0;j<2;j++){
			p = ls[j];
			q = ls[1-j];
			nk = nd[j];
		       	nl = nd[1-j];
			for(size_t i=0;i<p->size();i++){
				// find whether block_id in another list
				size_t block_id = (*p)[i];
				it_1 = find( (*q).begin(), (*q).end(), block_id);
				// 2 nodes in the same block
				if(it_1!=(*q).end() && nk->isS()!=X && 
				  !nk->is_ground() && (nk->nbr[TOP]!=NULL && 
				     nk->nbr[TOP]->type!=INDUCTANCE)){
				   //clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
					Block &block = block_info[block_id];
					size_t k1 = nk->id_in_block[i];
					nr = nl->nbr[TOP]->ab[0]->rep;
					if(nr->isS()!=Y)
					   nr = nl->nbr[TOP]->ab[1]->rep;
					//clog<<"nr: "<<*nr<<endl;
					block.bp[k1] += G *(nr->value);
					//clog<<"k1, G, bp: "<<k1<<" "<<G<<" "<<block.bp[k1]<<endl;
				}
			}
		}
	}
}

void Circuit::make_A_symmetric_block_tr(Tran &tran){
	NetList::iterator it;
	vector<size_t> *p, *q;
	Node *nk, *nl;

      for(int type = 0;type <NUM_NET_TYPE;type++){  
	 NetList & ns = net_set[type];
	 if(!(type == RESISTOR || type ==INDUCTANCE))
	    continue;
	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
	
		Net *net = *it;
		Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};	
		//clog<<"symmetric net: "<<*net<<endl;
		if(!(nd[0]->isS()==Y || nd[1]->isS()==Y)) continue;
		double G=0;
		if(net->type == RESISTOR)
		  G = 1./net->value;
                else if(net->type ==INDUCTANCE)
                   G = tran.step_t/(2*net->value);
		//clog<<"symmetric net again: "<<*net<<endl;

		vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
		vector<size_t >::const_iterator it_1;
		for(size_t j=0;j<2;j++){
			p = ls[j];
			q = ls[1-j];
			nk = nd[j]->rep;
		       	nl = nd[1-j]->rep;
			for(size_t i=0;i<p->size();i++){
				// find whether block_id in another list
				size_t block_id = (*p)[i];
				it_1 = find( (*q).begin(), (*q).end(), block_id);
				// 2 nodes in the same block
				if(it_1!=(*q).end() && nk->isS()==X && 
				  !nk->is_ground()&& (nk->nbr[TOP]!=NULL && 
				     nk->nbr[TOP]->type==INDUCTANCE)){ 
                                   //clog<<"nk, nl: "<<*nk<<" "<<*nl<<endl;
                                   Block &block = block_info[block_id];
                                   size_t k1 = nk->id_in_block[i];
                                   if(!nl->is_ground()){
                                      block.b_trp[k1] += G *(nl->value);
                                      //clog<<"nl->value, G, bp: "<<nl->value<<" "<<G<<" "<<block.b_trp[k1]<<endl;
                                   }
                                }
			}
		}
	}
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

void Circuit::assign_task_tree(int N_proc, int N){
   size_t num = N / N_proc;
   int res = N % N_proc;
   size_t  tot = 0;
   size_t num_temp = 0;
   clog<<"num, res: "<<num<<" "<<res<<endl;
   for(int i=0;i<N_proc;i++){
      my_id =i;
      num_temp = num;
      if(i < res)
         num_temp++;
      start_tree[i] = tot;
      tot += num_temp;
      end_tree[i]= tot;
      // clog<<"proc, N, start, end: "<<i<<" "<<N<<" "<<start_tree[i]<<" "<<end_tree[i]<<endl;
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
bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2){
  return (nd_1->value < nd_2->value); 
}
/*void Circuit::update_node_set_bx(){
   int id = 0;
   //clog<<"len_node_set_b. "<<pg.node_set_b.size();
   for(int i=0;i<pg.node_set_b.size();i++){
      id = pg.node_set_b[i];
      pg.node_set_b[i] = id_map[id];
      //clog<<"b_old, b_new: "<<id<<" "<<id_map[id]<<endl;
   }
   //clog<<endl;
   //clog<<"len_node_set_x. "<<pg.node_set_x.size()<<endl;
   for(int i=0;i<pg.node_set_x.size();i++){
      id = pg.node_set_x[i];
      pg.node_set_x[i] = id_map[id];
      //clog<<"x_old, x_new: "<<id<<" "<<id_map[id]<<endl;
   }
}*/


// each thread's solve function
void *Circuit::ptr_solve_FFS(double *X){
   clock_t t1, t2;
   double *Lx;
   int *Li, *Lp, *Lnz;
   int i=0, j=0, p, lnz, pend;
   double y;

   int n = L->n;
   bool *done;
   done = new bool [L->n];
   // assign xp[i] = bnewp[i];
   //start_ptr_assign_1();
   for(i=0;i<n;i++){
	done[i] = false;
	X[i] = bnewp[i];
   }

   Lp = static_cast<int *>(L->p);  
   Lx = static_cast<double*> (L->x);
   Li = static_cast<int*>(L->i) ;
   Lnz = static_cast<int *>(L->nz);

  t1 = clock();
  for(i=0; i<tree.size();i++){
  	j = tree[i]->value;
	// copy diagonal element
	
	//****** solve_col_FFS_pr function
	p = Lp [j] ;
	lnz = Lnz[j] ;
	pend = p + lnz ;

	y = X[j] ;
	if(L->is_ll == true){
		X[j] /= Lx [p] ;
	}

	//cout<<"j, X: "<<j<<" "<<X[j]<<endl;
	for ( ++p ; p < pend ; p++){	
		X[Li [p]] -= Lx [p] * y ;
	}
	done[j] = true;
  }
  t2 = clock();
  //clog<<"solve: "<<tree.size()<<" cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;

  t1 = clock();
  for(j=0;j<n;j++){
  	if(done[j] == true) continue;
	p = Lp [j] ;
	lnz = Lnz[j] ;
	pend = p + lnz ;

	y = X[j] ;
	if(L->is_ll == true){
		X[j] /= Lx [p] ;
	}

	for ( ++p ; p < pend ; p++){	
		X[Li [p]] -= Lx [p] * y ;
	}
	//done[j] = true;
  }
  //for(j=0;j<n;j++)
	//cout<<"FFS sol: "<<j<<" "<<X[j]<<endl;
  t2 = clock();
  //clog<<"FFS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;


  t1 = clock();
  // FBS solve
  for(j=n-1;j>=0; j--){
	  if(done[j] == true) continue;
	  /* get the start, end, and length of column j */
	  p = Lp [j] ;
	  lnz = Lnz [j] ;
	  pend = p + lnz ;
	  double y = X [j] ;
	  double d = Lx [p] ;
	  if(L->is_ll == false)
		  X[j] /= d ;
	  for (++p ; p < pend ; p++)
	  {
		  X[j] -= Lx [p] * X [Li [p]] ;
	  }
	  if(L->is_ll == true)
		  X [j] /=  d ;
  }

  for(i=0; i<tree.size();i++){
  	j = tree[i]->value;
	// copy diagonal element
	
	//****** solve_col_FFS_pr function
	p = Lp [j] ;
	lnz = Lnz[j] ;
	pend = p + lnz ;

	y = X[j] ;
	if(L->is_ll == true){
		X[j] /= Lx [p] ;
	}

	for ( ++p ; p < pend ; p++){	
		X[Li [p]] -= Lx [p] * y ;
	}
	done[j] = true;
  }
  t2 = clock();
  //clog<<"FBS cost: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
  delete[] done;
  return 0;
}


void Circuit::solve_eq(double *X){
   int p, q, r, lnz, pend;
   int j, k, n = L->n ;
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

         double y = X [j] ;
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
   int j, k, n = L->n ;
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
         double y = X [j] ;
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

// solve rest columns with FFS and FBS
void Circuit::solve_single_col(cholmod_factor*L, double *X, double *Lx,
	int *Li, int *Lp, int *Lnz){
#if 0    
int i, j, n = L->n ;   
    //********* for the rest of nodes, solve in sequential
    Node_G *nd, *p;
    int col = 0;
    // solve the rest of node in FFS
    for(i = base_level[level_tr];i<tree.size();i++){	
	    // loop all the head nodes in level i
	    nd = tree[i];
	    p = nd;
	    // solve the list connected by head node
	    while(1){
		    col = p->value;
		    //clog<<"solve col: "<<col<<endl;
		    solve_col_FFS(L, X, col, Lx, Li, Lp, Lnz);
		    if(p->children == NULL) break;
		    if(p->children->level == p->level)
		    	p = p->children;
		    else break;
	    }
    }
    //for(i=0;i<n;i++)
	//cout<<"i, x after ffs: "<<i<<" "<<X[i]<<endl;
    // then solve these node in FBS
    vector<Node_G *> children;
    for(i = tree.size()-1; i >= base_level[level_tr]; i--){	
	    // loop all the head nodes in level i
	    nd = tree[i];
	    p = nd;
	    children.push_back(nd);
	    while(1){
		if(p->children == NULL) break;
		if(p->children->level == p->level){
			p = p->children;
			children.push_back(p);
		}
		else break;	
	    }

	    // solve the list connected by head node
	    for(j=children.size()-1; j>=0; j--){
		    p = children[j];
		    col = p->value;
		    solve_col_FBS(L, X, col, Lx, Li, Lp, Lnz);
	    }
	    children.clear();
    }
#endif
}
//#endif

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

   double y = X [j] ;
   double d = Lx [p] ;
   if(L->is_ll == false)
	   X[j] /= d ;
   for (++p ; p < pend ; p++)
	   X[j] -= Lx [p] * X [Li [p]] ;
   if(L->is_ll == true)
	   X [j] /=  d ;
}

void Circuit::build_etree(cholmod_factor *L, vector<Node_G*> &etree){
   double *Lx;
   int *Li, *Lp, *Lnz;
   int p, q, r, lnz, pend;
   Lp = static_cast<int *>(L->p);
   Lx = static_cast<double*> (L->x);
   Li = static_cast<int*>(L->i) ;
   Lnz = static_cast<int *>(L->nz);
   int j, k, n = L->n ;

   etree.clear();
   // first produce all nodes
   Node_G *nd;
   for(j=0;j<n;j++){
	nd = new Node_G(j);
	etree.push_back(nd);
  }

  for(j=0;j<n;j++){
  	p = Lp[j];
	lnz = Lnz[j];
	pend = p + lnz;
	nd = etree[Li[p]];

	if(++p < pend){
		nd->children = etree[Li[p]];
		etree[Li[p]]->parent.push_back(nd);
	} 
  }
  // assign super nodes level
  build_tree(etree); 
}

// locate root node, and find max_depth
void Circuit::build_tree(vector<Node_G*> &etree){
	Node_G *nd;
	int j;
	int n = replist.size();

	// first positive scan
	for(j=0;j<n;j++){
		nd = etree[j];
		if(nd->parent.size()!=0) continue;
		tree.push_back(nd);	
	}	
	//cout<<"leaf node size / n: "<<tree.size()<<" "<<n<<endl;	
}

// update level info between 2 nodes: nd and its parent
int Circuit::find_level(Node_G *nd){
	if(nd->children == NULL) {
		//tree.push_back(nd);
		return 1;
	}
	int level=0;
	//cout<<"parent_size: "<<nd->children->parent.size()<<endl;
	if(nd->children->parent.size()>1){
		level = nd->level +1;
		//tree.push_back(nd);
	}
	else if(nd->children->parent.size()==1){
		level = nd->level;
	}
	int diff = level - nd->children->level;
	//cout<<"nd, nd_l: "<<*nd<<endl;;
	//cout<<"parent, parent_l: "<<*nd->children; 
	// update parent level
	if(diff >0)
		nd->children->level += diff;	
	return 0;
}

// update level info between 2 nodes: nd and its parent
void Circuit::find_level_inv(Node_G *p, vector<Node_G*> &tree){
	if(p->parent.size() == 0){
		tree.push_back(p);
		return;
	}
	//int level=0;

	if(p->parent.size()>1){
		//level = p->level -1;
		tree.push_back(p);
	}
#if 0
	else if(p->parent.size()==1){
		level = p->level;
	}
#endif
	//clog<<"level: "<<level<<endl;
	//clog<<"parent size: "<<p->parent.size()<<endl;
#if 0	
	for(int j=0;j<p->parent.size();j++){
		p->parent[j]->level = level;
		//clog<<"parent node: "<<*p->parent[j];
	}
#endif
	for(int j=0;j<p->parent.size();j++){
		Node_G *q;
		q = p->parent[j];
		//p = p->parent[j];
		find_level_inv(q, tree);
	}
}

bool compare_s_level(const Node_G *nd_1, const Node_G *nd_2){
  return (nd_1->level < nd_2->level); 
}

#if 0
// FFS solve
void Circuit::solve_col_FFS_pr(cholmod_factor *L, double *X, int &j,
	double *Lx, int *Li, int *Lp, int *Lnz){
   int p, lnz, pend;
   
   p = Lp [j] ;
   lnz = Lnz[j] ;
   pend = p + lnz ;

   double y = X [j] ;
   if(L->is_ll == true)
	   X[j] /= Lx [p] ;
   for (p++ ; p < pend ; p++)
   	#pragma omp atomic
	   X [Li [p]] -= Lx [p] * y ;
}

// FBS solve
void Circuit::solve_col_FBS_pr(cholmod_factor *L, double *X, int &j,
	double *Lx, int *Li, int *Lp, int *Lnz){
   int p, lnz, pend;
   
   p = Lp [j] ;
   lnz = Lnz [j] ;
   pend = p + lnz ;

   double y = X [j] ;
   double d = Lx [p] ;
   if(L->is_ll == false)
	   X[j] /= d ;
   for (p++ ; p < pend ; p++)
	#pragma omp atomic
	   X[j] -= Lx [p] * X [Li [p]] ;
   if(L->is_ll == true)
	   X [j] /=  d ;
}
#endif

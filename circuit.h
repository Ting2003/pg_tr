// ----------------------------------------------------------------//
// Filename : circuit.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Circuit class
// used to construct the circuit network
// ----------------------------------------------------------------//
// - Ting Yu - Tue Feb 8 5:45 pm 2011
//   * added the ostream<< func in .h file
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __CIRCUIT_H__
#define __CIRCUIT_H__

#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <map>
#include <cmath>
#include "global.h"
#include "node.h"
#include "net.h"
#include "vec.h"
#include "triplet.h"
#include "block.h"
#include "transient.h"
#include "cholmod.h"
#include "etree.h"
#include <algorithm>

using namespace std;
using namespace std::tr1;

typedef vector<double> DoubleVector;
typedef vector<Net *> NetPtrVector;
typedef vector<Node *> NodePtrVector;
typedef NetPtrVector NetList;

// functor of translating Node * to void *
namespace std{ 
	namespace tr1{
		template<> struct hash< Node * >{
			size_t operator()( const Node * x ) const {
				return hash< const char* >()( (char*) x );
			}
		};
	}
}

class Circuit{
public:
	Circuit(string name="");
	~Circuit();
	void check_sys() const;
	friend class Block;
	// can be written as inline to speed up
	Node * get_node(string name);
	Net * get_net(string name);
	string get_name() const;

	static size_t get_total_num_layer();

	// add a node into nodelist
	bool add_node(Node * nd);

	// add a net into netset
	bool add_net(Net * net);

	bool has_node(string name) const;
	bool has_net(string name) const;

	// sort nodes according to predefined order
	void sort_nodes();
        
        void make_A_symmetric(double *bp);
	void make_A_symmetric_tr(double *b, double *x, Tran &tran);
	void make_A_symmetric_block();
      void make_A_symmetric_block_tr(Tran &tran);

      // solve for node voltage
	void solve(Tran &tran);
	
	void set_blocklist(Node * nd);

	static void set_parameters(double, double, double, size_t, int);
	static void get_parameters(double&, double&, double&, size_t&, int&);

	friend ostream & operator << (ostream & os, const Circuit & ckt);
	friend class Parser;

	// C style output
	void print();

private:
	// member functions
	void solve_LU(Tran &tran);
	void solve_LU_core(Tran &tran);

	bool solve_IT(Tran &tran);
	void solve_block_LU();

	bool solve_pcg();
	//bool solve_block_pcg();

	
	// initialize things before solve_iteration
	void solve_init(Tran &tran);
	void count_merge_nodes();

	// updates nodes value in each iteration
	double solve_iteration();
	bool solve_IT_dc(Matrix *A);
	bool solve_IT_tr(Matrix *A, Tran &tran);

	void block_init();
	void update_block_geometry();

	// methods of stamping the matrix
	void stamp_by_set(Matrix & A, double * b);
	void stamp_resistor(Matrix & A, Net * net);
	void stamp_current(double * b, Net * net);
	void stamp_VDD(Matrix & A, double *b, Net * net);
	void stamp_VDD_tr(double *b, Net * net);
	void stamp_inductance_dc(Matrix & A, double *b, Net * net);
	void stamp_capacitance_dc(Matrix & A, Net * net);

	void stamp_by_set_tr(Matrix & A, double *b, Tran &tran);
	void stamp_resistor_tr(Matrix & A, Net * net);
	void current_tr(Net *net, double &time);
	
	void stamp_current_tr_1(double *bp, double *b, Tran &tran, double &time);
	void stamp_current_tr_net_1(double *bp, double *b, Net *net, double &time);

	void stamp_current_tr(double *b, Tran &tran, double &time);
	void stamp_current_tr_net(double *b, Net * net, double &time);
	void stamp_capacitance_tr(Matrix & A, Net * net, Tran &tran);
	void stamp_inductance_tr(Matrix & A, Net * net, Tran &tran);
	void modify_rhs_tr(double *b, double *xp, Tran &tran, int &iter);
	void modify_rhs_c_tr(Net *net, double *rhs, double *xp, Tran &tran, int &iter);
	void modify_rhs_l_tr(Net *net, double *rhs, double *xp, Tran &tran, int &iter);
	void release_tr_nodes(Tran &tran);
	void link_tr_nodes(Tran &tran);
	void save_tr_nodes(Tran &tran, double *x);
	void print_tr_nodes(Tran &tran);

	void stamp_block_matrix(Matrix *A);
	void stamp_block_matrix_tr(Matrix *A, Tran &tran);

	void stamp_block_resistor(Net *net, Matrix * A);
	void stamp_block_current(Net * net);
	void stamp_block_VDD(Net * net, Matrix * A);
	void stamp_block_VDD_tr(Net * net);
	void stamp_block_inductance_dc(Net * net, Matrix * A);
	void update_block_rhs(Block & block, int dir);

	// transient stamp
	void stamp_block_resistor_tr(Net *net, Matrix * A);
	void stamp_block_current_tr(double &time);
	void stamp_block_current_tr_net(Net * net, double &time);
	void modify_block_rhs_tr(Tran &tran);

	void modify_block_rhs_c_tr(Net *net, Tran &tran);
	void modify_block_rhs_l_tr(Net *net, Tran &tran);

	void stamp_block_capacitance_dc(Net * net, Matrix * A);
	void stamp_block_resistor_cl_tr(Matrix *A, Net * net, 
			double &Geq);
	double solve_iteration_tr(Tran &tran, double &time);
	double solve_iteration_tr_step();
	//  ******* method for PCG method  ********
	// solve circuit with preconditioned pcg method
	void copy_node_voltages(double *x, size_t &size, bool from=true);
	void copy_node_voltages_block(bool from=true);

	// after solving, copy node voltage from replist to nodes
	void get_voltages_from_LU_sol(double *x);
	void get_voltages_from_block_LU_sol();

        void set_len_per_block();
	void find_block_size ();
	void block_boundary_insert_net(Net * net);
	void find_block_base();

	void partition_circuit();
	double modify_voltage(Block & block);
	double modify_voltage(Block & block, double *x_old);

	void node_voltage_init();
	void solve_one_block(size_t block_id);

	void select_omega();

	void set_type(CIRCUIT_TYPE type){circuit_type = type;};
        // ************* functions and members for thread **********

        double *temp;
	void assign_task(int N_proc, int N);

        void start_thread();
        static void* call_thread_task(void*arg);
        void *thread_task();

        void start_ptr_memset();
        static void *call_ptr_memset(void *arg);
        void *ptr_memset();

        void start_ptr_assign();
        static void *call_ptr_assign(void *arg);
        void *ptr_assign();

	void start_ptr_assign_rid();
        static void *call_ptr_assign_rid(void *arg);
        void *ptr_assign_rid();

	void start_ptr_solve();
        static void *call_ptr_solve(void *arg);
        void *ptr_solve();

	void start_ptr_assign_bp();
        static void *call_ptr_assign_bp(void *arg);
        void *ptr_assign_bp();

	void start_ptr_assign_xp();
        static void *call_ptr_assign_xp(void *arg);
        void *ptr_assign_xp();

	void start_ptr_assign_xp_b();
        static void *call_ptr_assign_xp_b(void *arg);
        void *ptr_assign_xp_b();

	void start_ptr_assign_1();
        static void *call_ptr_assign_1(void *arg);
        void *ptr_assign_1();


        pthread_t tid[NTHREADS];
        size_t start[NTHREADS];
        size_t end[NTHREADS];
        int my_id;

        int *id_map;
        cholmod_factor *L;
        cholmod_common c, *cm;
        cholmod_dense *b, *x, *bnew;
        double *bp, *xp;
        double *bnewp;
        void solve_eq(cholmod_factor *L, double *X);
        //void solve_eq_sp(cholmod_factor *L, double *X);
	void solve_eq_pr(cholmod_factor *L, double *X);
        void solve_col_FFS(cholmod_factor *L, double *X, int &j,
	double *Lx, int *Li, int *Lp, int *Lnz);
	
	void solve_col_FBS(cholmod_factor *L, double *X, int &j,
	double *Lx, int *Li, int *Lp, int *Lnz);
	// ************ elimination tree  ************
	vector<Node_G*> etree;
	vector<Node_G *> tree;
	int *n_level;
	int *base_level;
	void build_etree(cholmod_factor *L, vector<Node_G*> &etree);
	int find_level(Node_G *nd);
	void find_level_inv(Node_G *nd, vector<Node_G*> &tree);
	void build_tree(vector<Node_G*> &etree);
#if 0	
	// ****************functions and members for sp_v methods ***
        Path_Graph pg;
        int *path_b, *path_x;
        int len_path_b, len_path_x;
        
	void update_node_set_bx();
        void parse_path_table();
        void build_path_graph(cholmod_factor *L);
        void build_FFS_path(cholmod_factor *L);
        void build_FBS_path(cholmod_factor *L);
        void set_up_path_table(cholmod_factor *L);
        void find_path(vector<size_t>&node_set, List_G &path);
        
#endif
	// ************** member variables *******************
	NodePtrVector nodelist;		// a set of nodes
	NodePtrVector replist;		// a set of representative nodes
	NetPtrVector net_set[NUM_NET_TYPE];// should be the same as size of NET_TYPE
	// defines the net direction in layers
	static vector<LAYER_DIR> layer_dir;
	vector<int> layers;
	
	// mapping from name to Node object pointer
	unordered_map<string, Node*> map_node;

	// mapping from Net name to object pointer
	// unordered_map<string, Net*> map_net;

	// mapping from Node pointer to their index in nodelist
	unordered_map<Node *, size_t> node_id;
	//unordered_map<Node *, size_t> rep_id;

	// circuit name
	string name;

	// blocks
	BlockInfo block_info;
	size_t x_min, y_min, x_max, y_max;

	// control variables
	static double EPSILON;
	static double OMEGA;
	static double OVERLAP_RATIO;
	static size_t MAX_BLOCK_NODES;
	static int MODE; // 0 = IT, 1 = LU

	CIRCUIT_TYPE circuit_type;

	double VDD;
	size_t num_blocks;
};

inline size_t Circuit::get_total_num_layer(){return layer_dir.size();}

// adds a node into nodelist
inline bool Circuit::add_node(Node * node){
	nodelist.push_back(node);
	map_node[node->name] = node;
	return true;
}

// adds a net into netset
inline bool Circuit::add_net(Net * net){
	//map_net[net->name] = net;
	net_set[net->type].push_back(net);
	return true;
}

// fina a node by name
inline bool Circuit::has_node(string name) const{
	if( map_node.find(name) != map_node.end() ) return true;
	return false;
}

// get a node by name
inline Node * Circuit::get_node(string name){
	unordered_map<string, Node*>::const_iterator it = map_node.find(name);
	if( it != map_node.end() ) return it->second;
	else return NULL;
}

/*
// find a net by name
inline bool Circuit::has_net(string name) const{
	if( map_net.find(name) != map_net.end() ) return true;
	return false;
}


// get a net by name
inline Net * Circuit::get_net(string name){return map_net[name];}
*/

bool compare_Node_G(const Node_G *nd_1, const Node_G *nd_2);
bool compare_s_level(const Node_G *nd_1, const Node_G *nd_2);

bool compare_node_ptr(const Node *a, const Node *b);
ostream & operator << (ostream & os, const NodePtrVector & nodelist);
ostream & operator << (ostream & os, const NetPtrVector & nets);
//ostream & operator << (ostream & os, const vector<Block > & block_info);
#endif

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
#include "transient.h"
#include "cholmod.h"
#include <algorithm>
#include "sp_graph_table.h"
#include "sp_node.h"

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
	//void make_A_symmetric_block();
      //void make_A_symmetric_block_tr(Tran &tran);

      // solve for node voltage
	void solve(Tran &tran);
	
	//void set_blocklist(Node * nd);
	friend ostream & operator << (ostream & os, const Circuit & ckt);
	friend class Parser;

	// C style output
	void print();

private:
	// member functions
	void solve_LU(Tran &tran);
	void solve_LU_core(Tran &tran);
	
	// initialize things before solve_iteration
	void solve_init();
	void count_merge_nodes();

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
	
	void stamp_current_tr_1(double *bp, double *b, double &time);
	void stamp_current_tr_net_1(double *bp, double *b, Net *net, double &time);

	void stamp_current_tr(double *b, double &time);
	void stamp_current_tr_net(double *b, Net * net, double &time);
	void stamp_capacitance_tr(Matrix & A, Net * net, Tran &tran);
	void stamp_inductance_tr(Matrix & A, Net * net, Tran &tran);
	void modify_rhs_tr_0(double *b, double *xp);

	void modify_rhs_tr(double *b, double *xp);
	void set_eq_induc(Tran &tran);
	void set_eq_capac(Tran &tran);
	void modify_rhs_c_tr_0(Net *net, double *rhs, double *xp);
	void modify_rhs_l_tr_0(Net *net, double *rhs, double *xp);

	void modify_rhs_c_tr(Net *net, double *rhs, double *xp);
	void modify_rhs_l_tr(Net *net, double *rhs, double *xp);
	void release_tr_nodes(Tran &tran);
	void release_ckt_nodes(Tran &tran);
	void print_ckt_nodes(Tran &tran);
	void save_ckt_nodes_to_tr(Tran &tran);
	void link_tr_nodes(Tran &tran);
	void link_ckt_nodes(Tran &tran);
	void save_tr_nodes(Tran &tran, double *x);
	void save_ckt_nodes(Tran &tran, double *x);

	void print_tr_nodes(Tran &tran);

	void copy_node_voltages(double *x, size_t &size, bool from=true);

	// after solving, copy node voltage from replist to nodes
	void get_voltages_from_LU_sol(double *x);
	void select_omega();

	void set_type(CIRCUIT_TYPE type){circuit_type = type;};
        // ************* functions and members for thread **********

        double *temp;	
        int *id_map;
        cholmod_factor *L;
	double *Lx;
	int *Li, *Lp, *Lnz;
        cholmod_common c, *cm;

        cholmod_dense *b, *x, *bnew;
        double *bp, *xp;
        double *bnewp;

	void solve_eq(double *X);
	void solve_eq_sp(double*X);
	// set s_col_FFS and FBS
	void solve_eq_set();
	int* s_col_FFS;
	int* s_col_FBS;
	// ********* sparse vectors ******
	Path_Graph pg;
	int *path_b, *path_x;
	int len_path_b, len_path_x;
	int flag_ck;
	 void find_super();
	 void update_node_set_bx();                               
         void parse_path_table();
         void build_path_graph();                
         void build_FFS_path();
         void build_FBS_path();                  
         void set_up_path_table();               
         void find_path(vector<size_t>&node_set, List_G &path);
       
	vector<Node_TR_PRINT> ckt_nodes;
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
	//BlockInfo block_info;
	//size_t x_min, y_min, x_max, y_max;

	// control variables

	CIRCUIT_TYPE circuit_type;

	double VDD;
	//size_t num_blocks;
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


ostream & operator << (ostream & os, const NodePtrVector & nodelist);
ostream & operator << (ostream & os, const NetPtrVector & nets);
//ostream & operator << (ostream & os, const vector<Block > & block_info);
#endif

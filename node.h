#ifndef __NODE_H__
#define __NODE_H__

#include <string>
#include <algorithm>
#include "global.h"
#include "point.h"
#include "net.h"
#include <vector>
#include <iostream>
//#include "block.h"
using namespace std;

class Net;
class Circuit;
// a Node in the network
class Node{
public:
	// member functions
	Node();
	Node(string name, Point _pt, int flag=-1, double v=0.0);
	void set_nbr(DIRECTION dir, Net * name);
	Net * get_nbr(DIRECTION dir) const;

	// Ting: get_nbr_node
	Node * get_nbr_node(Node *node, DIRECTION dir) const;
	Node get_nbr_node(Node node, DIRECTION dir) const;

	int get_layer() const;

	int isS() const;
	bool is_ground() const;

	double get_value() const;
	void set_value(double v);

	friend ostream & operator << (ostream & os, const Node & node);
	friend class Circuit;
	friend class Parser;

	////////////////////////////////////////////////////////////
	// member variables
	string name;		// node name
	Point pt;		// coordinate
	// only 2 possible cases:
	// {TOP, BOTTOM, EAST, WEST}
	// {TOP, BOTTOM, NORTH, SOUTH}
	Net * nbr[6];		// neighboring nets

	size_t rid;		// id in rep_list

private:
	double value;		// voltage
	// flag = 1 --> X
	// flag = 2 --> Y
	// flag = 3 --> Z
	int flag;		// mark if the node is an XYZ node
	Node * rep;		// representative, if somewhere is short-circuit
	//vector<size_t> blocklist;	// belongs to which block
	//vector<size_t> id_in_block;	// local index inside block	
};      	

inline int Node::isS() const{return flag;}

//inline bool Node::is_ground() const{return name == "0";}
// use a tricky way to speed up
inline bool Node::is_ground() const{return pt.x<0;}

inline int Node::get_layer() const{ return pt.z; }

inline double Node::get_value() const{return value;}

inline void Node::set_value(double v){value = v;}

inline void Node::set_nbr(DIRECTION dir, Net * net){ nbr[dir] = net; }

#endif

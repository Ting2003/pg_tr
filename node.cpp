#include <iomanip>
#include "global.h"
#include "node.h"
using namespace std;

// empty constructor
Node::Node():name(""),pt(Point(-1,-1,-1)), rid(0),
	value(0.0), flag(-1), rep(NULL){
	for(int i=0;i<6;i++) this->nbr[i] = NULL;
}

Node::Node(string n, Point _pt, int x, double v): 
	name(n), pt(_pt), rid(0), 
	value(v), flag(x), rep(NULL) {
	for(int i=0;i<6;i++) this->nbr[i] = NULL;
}

// Ting
// get node's nbr from dir direction
// template<class T>
Node * Node::get_nbr_node(Node *node, DIRECTION dir) const{
	Node * node_nbr =node->nbr[dir]->ab[0];
	if(node_nbr->pt == node->pt){
		node_nbr = node->nbr[dir]->ab[1];
	}
	return node_nbr;	
}

// Ting
// get nodes' nbr from dir direction
// template<class T>
Node Node::get_nbr_node(Node node, DIRECTION dir) const{
	Node node_nbr = *node.nbr[dir]->ab[0];
	if(node_nbr.pt == node.pt){
		node_nbr = *node.nbr[dir]->ab[1];
	}
	return node_nbr;
}

ostream & operator << (ostream & os, const Node & node){
	//os<<node.name<<node.pt<<"="<<scientific<<node.value;
	//os<<" rep="<<node.rep->name;
	//Net * net = node.nbr[TOP];
	os<<setw(OUTPUT_WIDTH_STRING)<<node.name
//	  <<" -> "<<setw(OUTPUT_WIDTH_STRING)<<node.rep->name
	  //<<" top="<<(net==NULL?"NULL":net->name)
	  <<setw(OUTPUT_WIDTH_FLOAT)<<scientific<<node.value;
	return os;
}

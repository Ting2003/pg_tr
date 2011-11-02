#include <iostream>
#include <iomanip>
#include "node.h"
#include "net.h"
using namespace std;

// a, b are pointers to two nodes
//Net::Net(NET_TYPE t, string n, double v, Node * a, Node * b): 
//type(t), name(n), value(v){
Net::Net(NET_TYPE t, double v, Node * a, Node * b):
	type(t),value(v){
	ab[0]=a;
	ab[1]=b;
	V1=0; V2=0; TD=0; Tr=0; Tf=0; PW =0; Period =0; 
}

ostream & operator << (ostream & os, const Net & net){
	os//<<net.name
		<<"("
	    <<net.ab[0]->name<<","
	    <<net.ab[1]->name<<")="
 	    <<scientific<<net.value;
	return os;
}

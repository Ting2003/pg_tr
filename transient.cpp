#include "transient.h"
#include <stdlib.h>
#include <stdio.h>

Node_TR_PRINT::Node_TR_PRINT(){
	value.clear();
	node = NULL;
        flag = -1;
}

Node_TR_PRINT::~Node_TR_PRINT(){
	value.clear();
	node = NULL;
	flag = -1;
}

Tran::Tran(){
	step_t = 0;
	tot_t = 0;
	length = 0;
	isTran = 0;
}
Tran::~Tran(){
	nodes.clear();
}

// print out transient solutions
// time in ns, value in V
void Tran:: print_tr_nodes(){
	double time = 0;
	size_t j=0;
	for(size_t i=0;i<nodes.size();i++){
		time = 0;
		j=0;
		cout<<endl<<"Node: "<<nodes[i].name<<endl<<endl;
		while(time < tot_t){// && iter <1){
			printf(" %.3e %.6e\n", time, 
				nodes[i].value[j]);
			j++;
			time += step_t;
		}
		cout<<"END: "<<nodes[i].name<<endl;
	}	
}

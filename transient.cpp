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
	int iter = 0;
	for(size_t i=0;i<nodes.size();i++){
		time = 0;
		iter = 0; 
		j=0;
		cout<<"NODE: "<<nodes[i].name<<endl;
		while(time < tot_t*1e9){// && iter <1){
			printf("%.5e  %.5e\n", time, 
				nodes[i].value[j]);
			j++;
			iter ++;
			time += step_t*1e9;
		}
		cout<<"END: "<<nodes[i].name<<endl;
	}	
}

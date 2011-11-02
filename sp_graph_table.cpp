#include "sp_graph_table.h"

Path_Graph::~Path_Graph(){
   node_set_b.clear();
   node_set_x.clear();
   nodelist.clear();
}
ostream & operator << (ostream & os, const vector<int> &path){
   for(size_t i=0;i<path.size();i++)
      clog<<path[i]<<" ";
   clog<<endl;
   return os;
}

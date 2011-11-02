#ifndef _PATH_GRAPH_H
#define _PATH_GRAPH_H
#include "sp_global.h"
#include "sp_node.h"

// node in path graph
class Path_Graph{
public:
   vector<size_t> node_set_b;
   vector<size_t> node_set_x;
   vector<Node_G*> nodelist;
   List_G path_FFS;
   List_G path_FBS;

   ~Path_Graph();
   friend ostream & operator<<(ostream &os, const vector<int> &path); 
};

#endif

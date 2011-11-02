#ifndef _SP_NODE_H_
#define _SP_NODE_H_
#include "sp_global.h"

// node in path graph
class Node_G{
public:
   int flag; // flag 
   int value; // stores the id
   Node_G * next;

   Node_G();
   friend ostream & operator<<(ostream &os, const Node_G &node);
   bool is_eq(const Node_G *nb);
};

// list for Node_G, to build path_list
class List_G{
public:
   Node_G * first;
   Node_G * last;
   int size;

   List_G();
   List_G(Node_G *node);
   void destroy_list();
   void add_node(Node_G *node);
   Node_G * insert_node(Node_G *node, Node_G *nd);
   friend ostream & operator<<(ostream &os, const List_G *list);
   void assign_size();
   int get_size();

};

#endif

#ifndef _SP_NODE_H_
#define _SP_NODE_H_
#include "global.h"

// node in path graph
class Node_G{
public:
   int flag; // flag 
   int value; // stores the id
   Node_G * next;
   Node_G * previous;

   Node_G();
   friend ostream & operator<<(ostream &os, const Node_G &node);
   bool is_eq(const Node_G *nb);
};

// node in matrix list
class Node_A{
public:
   int i; // flag 
   int j;
   double x;
   Node_A * next;
   Node_A * previous;

   Node_A();
   Node_A(int row, int col, double value);
   friend ostream & operator<<(ostream &os, const Node_A &node);
   bool is_eq(const Node_A *nb);
};

// list for Node_G, to build path_list
class List_G{
public:
   Node_G * first;
   Node_G * last;
   int size;

   List_G();
   List_G(Node_G *node);
   void add_node(Node_G *node);
   void insert_node(Node_G *node);
   friend ostream & operator<<(ostream &os, const List_G *list);
   void assign_size();
   int get_size();

};

class List{
public:
   Node_A * first;
   Node_A * last;
   int size;

   List();
   List(Node_A *node);
   int get_size();
   void add_node(Node_A *node);
   friend ostream & operator<<(ostream &os, const List *list);

};

// node in path graph
class Table_Item{
public:
   int start; // flag 
   int end; // stores the id

   Table_Item();
};
#endif

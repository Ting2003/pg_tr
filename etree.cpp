#include "etree.h"
using namespace std;

Node_G::Node_G(){
   value = 0;
   n_child = 0;
   level = 0;
   parent=NULL;
}
ostream & operator << (ostream & os, const Node_G & node){
      os <<node.value<<" "<<node.n_child<<" "<<node.level<<endl;  
      return os;
}

bool Node_G::is_eq( const Node_G *nb){
   if(value == nb->value && n_child == nb->n_child && level == nb->level)
      return true;
   else
      return false; 
}

List_G::List_G(){
   first = NULL;
   last = NULL;
   size = 0;
}

List_G::List_G(Node_G *node){
   first = node;
   last = node;
}

void List_G::add_node(Node_G *node){
   if (first ==NULL){
      first = node;
      last = node;
   }
   else{

      last->parent = node;
      last = node;
   }
}

Node_G * List_G::insert_node(Node_G *node, Node_G *nd){
   // p points to previous node of nd
   Node_G *p;
   p = nd;
   
   while(nd->value <node ->value){
      if(nd->parent == NULL) break;
         if(p != nd ) p = p->parent;
         nd = nd->parent;
   }
   // first is the minimum
   if(nd->is_eq(first)){
    nd->parent = node;
    p = node; 
   }
   else{
      p->parent = node;
      node->parent = nd;
      p = node;
   }
   return p;
}

int List_G::get_size(){
   return size;
}

void List_G::destroy_list(){
   first = NULL;
   last = NULL;
   size = 0;
}

void List_G::assign_size(){
   Node_G *nd;
   nd = first;
   size =0;
   while(!(nd->is_eq(last))){
      size++;
      nd = nd->parent;
   }
   size++;
}

ostream & operator << (ostream & os, const List_G * list){
   Node_G *nd;
   nd = list->first;
   while(!(nd->is_eq(list->last))){
      os << *nd;
      nd = nd->parent;
   }
   os <<*(list->last)<<endl;
   return os;
}

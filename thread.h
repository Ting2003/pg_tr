#ifndef _THREAD_H_
#define _THREAD_H_
#include "circuit.h"
#include "globa.h"
#include "node.h"

class Thread{
public:
   int my_id;
   Circuit *cir;
   // shows the start and end id for task
   int *start;
   int *end;
   
   static void * call_thread_task(void*arg){
      return ((Thread*)arg)->thread_task();
   }

   void *thread_task();
   void start_thread();
}
#endif

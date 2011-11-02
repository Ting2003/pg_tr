#include "thread.h"

void *Thread::thread_task(){
   for(size_t i=start[my_id];i<end[my_id];i++){
      Node * p = cir->nodelist[i];
      cir->node_id[p] = i;
   }
   return 0;
}

void Thread::start_thread(){
   pthread_t tid[NTHREADS];
   start = new int [NTHREADS];
   end = new int [NTHREADS];
   int i=0;
   int result;
   assign_task(NTHREADS, cir.nodelist.size());
   for(i=0;i<NTHREADS;i++){
      my_id = i;
      //result = 
      pthread_create(&tid[i],NULL, Thread::call_thread_task, this);
   //if(result ==0)
      //pthread_detach(tid[i]);
      pthread_join(tid[i], NULL);
   }
   free(start);
   free(end);
}

void Thread::assign_task(int N_proc, int N){
   size_t num = N / N_proc;
   size_t  tot = 0;
   for(int i=0;i<N_prorc;i++){
      size_t num_temp = num;
      if(tot + num < N)
         num_temp++;
      start[i] = tot;
      tot += num_temp;
      end[i]= tot;
      clog<<"proc, N, start, end: "<<i<<" "<<N<<" "<<start[i]<<" "<<end[i]<<endl;
   }
}



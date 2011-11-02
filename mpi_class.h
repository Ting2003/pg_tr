#ifndef __MPI_CLASS_H_
#define __MPI_CLASS_H_

#include <iostream>

using namespace std;

class MPI_CLASS{
public:
	// member function
	MPI_CLASS();
	~MPI_CLASS();
	// member
	int NUM_NET_TYPE;
	int X_BLOCKS;
	int Y_BLOCKS;
	long x_max;
	long y_max;
	long x_min;
	long y_min;
	float overlap_ratio;
	double len_per_block_x;
	double len_per_block_y;
	double len_ovr_x;
	double len_ovr_y;

	int *start_task;
	int *end_task;
	int *tasks_n;

	float *geo;
	// block_geo after overlap enlargement
	float *block_geo;
	// block_geo origin, before overlap enlargement
	float *block_geo_origin;

	int block_size;

	// function
	void MPI_Assign_Task(int & num_procs);
	void set_geo_origin(MPI_CLASS &mpi_class);
};

#endif

#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <utility>
#include <limits>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "omp.h"
using std::pair;
using namespace std;
// important: each node has at most 4 connected nets
enum DIRECTION{WEST, EAST, SOUTH, NORTH, BOTTOM, TOP, UNDEFINED};
enum NET_TYPE{RESISTOR, CURRENT, VOLTAGE, CAPACITANCE, INDUCTANCE};
enum S_NODE{X, Y, Z};
const int NUM_NET_TYPE = 5;
enum LAYER_DIR{HR, VT, NA}; // NA means not available
enum RUN_MODE{IT,LU};
enum CIRCUIT_TYPE{WB, C4, UNKNOWN};

const size_t INFTY = std::numeric_limits<size_t>::max();
const int MAX_BUF = 512;
const int MAX_LAYER = 30;
const double VIRTUAL_RESISTER=1e-10;
const double FLOAT_ZERO = 1e-30;
const int OUTPUT_WIDTH_FLOAT = 15;
const int OUTPUT_WIDTH_INDEX = 3;
const int OUTPUT_WIDTH_STRING = 10;
// overlap_ratio should < 0.45, or will be divergence
//const size_t X_BLOCKS = 2; // # of blocks along x axis
//const size_t Y_BLOCKS = 1; // # of blocks along y axis
const int THRESHOLD = 1e6;
const int NTHREADS=16;
typedef pair<size_t,size_t> SizeTPair;
class Triplet;
typedef Triplet Matrix;
//typedef vector<Net *> NetPtrVector;
#endif

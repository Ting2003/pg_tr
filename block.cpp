// ----------------------------------------------------------------//
// Filename : block.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// implementation of block class
// ----------------------------------------------------------------//
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added opeartor !=
//   * added this log

#include <cassert>
#include "block.h"
#include "node.h"
#include "util.h"
#include "umfpack.h"

Block::Block(size_t _count):
        b_ck(NULL),
	b_new_ck(NULL),
	//bp(NULL),
	//bnewp(NULL),
	//xp(NULL),
	x_ck(NULL),	
	count(_count),
	nodes(NULL),
	lx(-1.0), ly(-1.0),
	ux(-1.0), uy(-1.0){}

Block::~Block(){
    delete [] nodes;
}

void Block::free_block_cholmod(cholmod_common *cm){
    cholmod_free_factor(&L, cm);
    cholmod_free_dense(&b_ck, cm);
    cholmod_free_dense(&b_tr, cm);
    cholmod_free_dense(&b_new_ck, cm);
    cholmod_free_dense(&x_ck, cm);
}

void Block::CK_decomp(Matrix & A, cholmod_common *cm){
	Algebra::CK_decomp(A, L, cm);
}

void Block::solve_CK(cholmod_common *cm){
	x_ck = cholmod_solve(CHOLMOD_A, L, b_new_ck, cm);
}

void Block::allocate_resource(cholmod_common *cm){
	if( count == 0 ) return;
	nodes = new Node *[count];

	b_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	x_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bp = static_cast<double*>(b_ck->x);
	xp = static_cast<double*>(x_ck->x);
	b_new_ck = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
	bnewp = static_cast<double*>(b_new_ck->x);    
	b_tr = cholmod_zeros(count, 1, CHOLMOD_REAL, cm);
        b_trp = static_cast<double *>(b_tr->x);
}

// return the relative position of this block to another block
DIRECTION Block::relative_to(const Block & block) const{
	if( bx == block.bx ){
		if( by - block.by == 1 )
			return NORTH;
		else if( block.by - by  == 1 )
			return SOUTH;
	}
	else if( by == block.by ){
		if( bx - block.bx == 1 )
			return EAST;
		else if ( block.bx - bx == 1 )
			return WEST;
	}
	//fprintf(stderr, "%ld %ld, %ld %ld\n",bx,by,block.bx,block.by);
	return UNDEFINED;	// add this to avoid compiler warning
}

// update rhs of each block with its boundary netlist
void Block::update_rhs(int flag_tr){
	size_t size = boundary_netlist.size();
	size_t k=0, l=0;
        if(flag_tr ==0){
           for(k=0;k<count;k++)
            bnewp[k] = bp[k];
        }
        else if(flag_tr ==1){
           for(k=0;k<count;k++)
            bnewp[k] = b_trp[k];
        }

	// for each net in this block
	for(size_t i=0;i<size;i++){
		Net * net = boundary_netlist[i];
		//clog<<"net: "<<*net<<endl;
		double G = 1.0/net->value;

		Node * a = net->ab[0]->rep;
		Node * b = net->ab[1]->rep;

		vector<size_t> &block_id_a = a->blocklist;
		vector<size_t> &block_id_b = b->blocklist;

		vector<size_t>::const_iterator it;

		// find out which end of the net is inside the block
		it = find(block_id_a.begin(), block_id_a.end(), bid);
		if(it != block_id_a.end()){
			k = a->id_in_block[it - block_id_a.begin()];
			//clog<<"a->id: "<<k<<endl;
			if(!a->is_ground() && a->isS()!=Y){
				bnewp[k] += G * b->value;
				//clog<<"G, b->value, b_new: "<<G<<" "<<b->value<<" "<<b_new[k]<<endl;
			}
		}
		else {
			it = find(block_id_b.begin(),
					block_id_b.end(), bid);
			if(it !=block_id_b.end() ){
				l = b->id_in_block[it - block_id_b.begin()];
				//clog<<"b->id: "<<l<<endl;
				if(!b->is_ground() && b->isS()!=Y){ 
					bnewp[l] += G *a->value;
					//clog<<"G, a->value, b_new: "<<G<<" "<<a->value<<" "<<b_new[k]<<endl;
				}
			}
		}
	} // end of for i
}

/////////////////////////////////////////////////////////////////
// methods for BlockInfo

Block * BlockInfo::get_block_neighbor(const Block & b, DIRECTION dir) {
	switch(dir){
	case WEST:
		if( b.bx >= 1 )
			return &blocks[b.bid-1];
		break;
	case EAST:
		if( b.bx + 1 < X_BLOCKS )
			return &blocks[b.bid+1];
		break;
	case SOUTH:
		if( b.by >= 1 )
			return &blocks[b.bid - X_BLOCKS];
		break;
	case NORTH:
		if( b.by + 1 < Y_BLOCKS)
			return &blocks[b.bid + X_BLOCKS];
		break;
	default:
		report_exit("Unkown direction\n");
		break;
	}
	// no neighbor
	return NULL;
}


void BlockInfo::update_block_geometry(){
	// compute the geometrical information for the blocks
	for(size_t y=0;y<Y_BLOCKS;y++){
		for(size_t x=0;x<X_BLOCKS;x++){
			size_t block_id = y * X_BLOCKS + x;
			Block & b = blocks[block_id];
			b.bx = x;
			b.by = y;
			b.bid = block_id;
			b.lx = x * len_per_block_x - len_ovr_x;
			b.ly = y * len_per_block_y - len_ovr_y;
			b.ux = (x+1) * len_per_block_x + len_ovr_x;
			b.uy = (y+1) * len_per_block_y + len_ovr_y;
		}
	}
}

void BlockInfo::set_len_per_block(size_t x_min, size_t x_max,
			          size_t y_min, size_t y_max,
				  double overlap_ratio){
	double x, y;
	//x = (double)(x_max-x_min) / X_BLOCKS;
	//y = (double)(y_max-y_min) / Y_BLOCKS;
	x = (double)(x_max-x_min+0.5) / X_BLOCKS;
	y = (double)(y_max-y_min+0.5) / Y_BLOCKS;
	//if( fzero(x) ) x = 1.0;
	//if( fzero(y) ) y = 1.0;
	len_per_block_x = x;
	len_per_block_y = y;
	len_ovr_x = x * overlap_ratio;
	len_ovr_y = y * overlap_ratio;
	clog<<"len_x, len_y: "<<x<<" / "<<y<<endl;
	clog<<"len_ovr_x, len_ovr_y: "<<len_ovr_x
		<<" / "<<len_ovr_y<<endl;
}

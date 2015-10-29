#ifndef H_CELL_FUNCTIONS
#define H_CELL_FUNCTIONS

#include "mesh.h"

int get_adjacent_index(int direction, int cur_y, int cur_x);
void compute_beta(cell_t *mesh, double beta_coef);

#endif /* H_CELL_FUNCTIONS */

#ifndef H_CELL_FUNCTIONS
#define H_CELL_FUNCTIONS

#include "mesh.h"

int get_adjacent_index(int direction, int cur_y, int cur_x);
void compute_beta(cell_t *mesh, int cur_y, int cur_x, double beta_coef);
void compute_A(cell_t *mesh, int cur_y, int cur_x);
void update_interior(cell_t *mesh, cell_t *mesh_old, int cur_y, int cur_x);
void update_boundary(cell_t *mesh, cell_t *mesh_old, int cur_y, int cur_x, int boundary_side);
void update_corner(cell_t *mesh, cell_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2);

#endif /* H_CELL_FUNCTIONS */

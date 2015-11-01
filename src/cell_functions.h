#ifndef H_CELL_FUNCTIONS
#define H_CELL_FUNCTIONS

#include "mesh.h"
#include "util.h"

int get_adjacent_index(mesh_t *mesh, int direction, int cur_y, int cur_x);
void compute_beta(mesh_t *mesh, int cur_y, int cur_x, double beta_coef);
void compute_A(mesh_t *mesh, int cur_y, int cur_x);
void update_interior(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x);
void update_boundary(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x, int boundary_side);
void update_corner(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2);

#endif /* H_CELL_FUNCTIONS */

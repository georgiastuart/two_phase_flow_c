#ifndef H_CELL_FUNCTIONS
#define H_CELL_FUNCTIONS

#include "util.h"

struct mesh;
typedef struct mesh mesh_t;

typedef struct cell
{
    /* Values that live at the center of each cell */
    double perm, pressure, source;

    /* Values that live along the edges */
    double flux[4], l[4], beta[4], robin[4], A[4];

    /* Final velocity data */
    double velocity_y, velocity_x;
} cell_t;

typedef struct cell_ops
{
    void (*cell_compute_beta)();
    void (*cell_compute_A)();
    void (*cell_update_interior)();
    void (*cell_update_boundary)();
    void (*cell_update_corner)();
} cell_ops_t;

const cell_ops_t cell_p_ops;

int get_adjacent_index(mesh_t *mesh, int direction, int cur_y, int cur_x);
void cell_p_compute_beta(mesh_t *mesh, int cur_y, int cur_x, double beta_coef);
void cell_p_compute_A(mesh_t *mesh, int cur_y, int cur_x);
void cell_p_update_interior(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x);
void cell_p_update_boundary(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side);
void cell_p_update_corner(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2);
void initialize_cell_ops(cell_ops_t *cell_ops, const char *mode);

#endif /* H_CELL_FUNCTIONS */

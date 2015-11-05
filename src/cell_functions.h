#ifndef H_CELL_FUNCTIONS
#define H_CELL_FUNCTIONS

#include "util.h"

struct mesh;
typedef struct mesh mesh_t;

typedef struct cell
{
    /* Values that live at the center of each cell */
    double perm, pressure, source, saturation, saturation_prev, diffusion, source_d;

    /* Values that live along the edges */
    double flux_p[4], l_p[4], beta[4], robin[4], A_p[4];
    double flux_d[4], l_d[4], beta_d[4], A_d[4];

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
    void (*cell_compute_robin)();
} cell_ops_t;

extern const cell_ops_t cell_press_ops, cell_diff_ops;

/* Gets adjacent mesh index */
int get_adjacent_index(mesh_t *mesh, int direction, int cur_y, int cur_x);

/* Beta, Robin, and A cell calculations */
void press_compute_beta(mesh_t *mesh, int cur_y, int cur_x, double beta_coef);
void diff_compute_beta(mesh_t *mesh, int cur_y, int cur_x, double beta_coef);
void press_compute_robin(mesh_t *mesh, int cur_y, int cur_x);
void diff_compute_robin(mesh_t *mesh, int cur_y, int cur_x);
void press_compute_A(mesh_t *mesh, int cur_y, int cur_x);
void diff_compute_A(mesh_t *mesh, int cur_y, int cur_x);

/* Diffusion coefficient calculation */
void diff_compute_diffusion(mesh_t *mesh, int cur_y, int cur_x);

/* Saturation or pressure updates */
void press_update_interior(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x);
void diff_update_interior(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x);
void press_update_boundary(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side);
void diff_update_boundary(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side);
void press_update_corner(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2);
void diff_update_corner(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side1, int boundary_side2);

/* For Diffusion test */
void diff_update_corner_dirichlet(mesh_t *mesh, mesh_t *mesh_old, int cur_y,
	int cur_x, int boundary_side1, int boundary_side2);
void diff_update_boundary_dirichlet(mesh_t *mesh, mesh_t *mesh_old, int cur_y, int cur_x,
                    int boundary_side);

#endif /* H_CELL_FUNCTIONS */

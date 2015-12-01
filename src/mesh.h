#ifndef H_MESH
#define H_MESH

#define MESH_INDEX(y, x) ((y + 1) * (mesh->dim.xdim + 2) + (x + 1))
#define MESH_INDEX_NO_PAD(y, x) (y * mesh->dim.xdim + x)
#define MESH_INDEX_INC_PAD(y, x) (y * (mesh->dim.xdim + 2) + x)

#include "util.h"
#include "cell_functions.h"
#include "mpi_util.h"

typedef struct global_mesh_params
{
    double porosity, visc_o, visc_w, sat_rel_o, sat_rel_w, eta, beta_coef;
} global_mesh_params_t;

typedef struct mesh
{
    cell_t *cell;
    dim_t dim;
    global_mesh_params_t global;
} mesh_t;

mesh_t* mesh_init_mesh(dim_t dim, double *perm, double *source, config_t *config);
void mesh_update(mesh_t *mesh, mesh_t *mesh_old, int block_type, const cell_ops_t *cell_ops);
int mesh_press_convergence_check(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff, int rank);
int mesh_diff_convergence_check(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff, int rank);
void mesh_press_impose_0_average(mesh_t *mesh, int rank);
void mesh_update_robin(mesh_t *mesh, const cell_ops_t *cell_ops);
void mesh_compute_velocity(mesh_t *mesh);
int mesh_pressure_iteration(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff,
        int block_type, int rank, send_vectors_t *send_vec, receive_vectors_t *rec_vec);
int mesh_diffusion_iteration(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff,
    int block_type, int rank, send_vectors_t *send_vec, receive_vectors_t *rec_vec);
void mesh_update_saturation_time(mesh_t *mesh, mesh_t *mesh_old);
void mesh_max_time_step(mesh_t *mesh, mesh_t *mesh_old);

/* For diffusion test only */
void setup_diffusion_test(mesh_t *mesh);

#endif /* H_MESH */

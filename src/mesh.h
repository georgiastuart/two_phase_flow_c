#ifndef H_MESH
#define H_MESH

#define MESH_INDEX(y, x) ((y + 1) * (mesh->dim.xdim + 2) + (x + 1))
#define MESH_INDEX_NO_PAD(y, x) (y * mesh->dim.xdim + x)
#define MESH_INDEX_INC_PAD(y, x) (y * (mesh->dim.xdim + 2) + x)
#define UP 0
#define RIGHT 1
#define DOWN 2
#define LEFT 3
#define ROBIN_MODE 0
#define SAT_MODE 1

#include "util.h"
#include "cell_functions.h"
#include "mpi_util.h"
#include "parameters.h"

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

typedef struct prod_well
{
    int x_pos, y_pos;
    double *oil_sat_recording;
} prod_well_t;

typedef struct production_wells
{
    int num_wells;
    prod_well_t *wells;
} production_wells_t;

mesh_t* mesh_init_mesh(dim_t dim, double *perm, double *source, double *sat, config_t *config);
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
void mesh_update_saturation_time(mesh_t *mesh);
void mesh_max_time_step(mesh_t *mesh, mesh_t *mesh_old);
int mesh_transport_iteration(mesh_t *mesh, mesh_t *mesh_old, int block_type, int rank,
            send_vectors_t *send_vec, receive_vectors_t *rec_vec);

/* Production Wells */
production_wells_t init_production_wells(mesh_t *mesh);
void record_production_wells(production_wells_t *wells, mesh_t *mesh, int time_step);

/* For tests only */
void setup_diffusion_test(mesh_t *mesh);
void setup_transport_test(mesh_t *mesh);

#endif /* H_MESH */

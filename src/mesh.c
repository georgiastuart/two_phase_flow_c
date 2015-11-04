#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "cell_functions.h"

/* Computes beta and A at all mesh points */
void mesh_compute_beta_A(mesh_t *mesh, const cell_ops_t *cell_ops)
{
    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cell_ops->cell_compute_beta(mesh, i, j, mesh->global.beta_coef);
            cell_ops->cell_compute_A(mesh, i, j);
        }
    }
}



/* Sets up the mesh with cell structs at each gridpoint */
mesh_t* mesh_init_mesh(dim_t dim, double *perm, double *source, config_t *config)
{
    int i, j;
    mesh_t *mesh;
    cell_t *cur_cell;

    /* Allocates memory for the mesh */
    mesh = malloc(sizeof(mesh_t));
    mesh->cell = malloc((dim.ydim + 2) * (dim.xdim + 2) * sizeof(cell_t));

    /* Sets the mesh dimensions */
    mesh->dim = dim;

    /* Sets global parameters */
    mesh->global.porosity = config->porosity;
    mesh->global.sat_rel_o = config->sat_rel_o;
    mesh->global.sat_rel_w = config->sat_rel_w;
    mesh->global.visc_o = config->visc_o;
    mesh->global.visc_w = config->visc_w;
    mesh->global.eta = config->eta;
    mesh->global.beta_coef = config->beta_coef;

    /* Sets cell permeability and sources */
    for (i = 0; i < (mesh->dim.ydim + 2); i++) {
        for (j = 0; j < (mesh->dim.xdim + 2); j++) {
            cur_cell = &mesh->cell[MESH_INDEX_INC_PAD(i, j)];
            cur_cell->perm = config->perm_scale * exp(config->perm_strength
                    * perm[MESH_INDEX_INC_PAD(i, j)]);
            cur_cell->source = source[MESH_INDEX_INC_PAD(i, j)];
            cur_cell->saturation = config->sat_rel_w + 0.01;
        }
    }

    /* Computes beta and A at all mesh points */
    mesh_compute_beta_A(mesh, &cell_p_ops);

    return mesh;
}

/* One iteration of the algorithm over all mesh points */
/* 0 - up, 1 - right,  2 - down, 3 - left */
void update_9(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Updates the corners */
    cell_ops->cell_update_corner(mesh, mesh_old, 0, 0, 0, 3);
    cell_ops->cell_update_corner(mesh, mesh_old, 0, mesh->dim.xdim - 1, 0, 1);
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, 0, 2, 3);
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, mesh->dim.xdim - 1, 2, 1);

    /* Updates the boundaries */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);
    }

    for (i = 1; i < (mesh->dim.xdim - 1); i++) {
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);
    }

    /* Updates the interior cells */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        for (j = 1; j < (mesh->dim.xdim - 1); j++) {
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* Iteration for a type 0 block */
void update_0(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, 0, 0, 0, 3);

    for (i = 1; i < mesh->dim.xdim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 1 block */
void update_1(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.xdim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 2 block */
void update_2(mesh_t *mesh, mesh_t *mesh_old, const const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, 0, mesh->dim.xdim - 1, 0, 1);

    for (i = 0; i < (mesh->dim.xdim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 3 block */
void update_3(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 4 block */
void update_4(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* Iteration for a type 5 block */
void update_5(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 6 block */
void update_6(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, 0, 2, 3);

    for (i = 1; i < mesh->dim.xdim; i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);

    for (i = 0; i < (mesh->dim.ydim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 7 block */
void update_7(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    for (i = 0; i < mesh->dim.xdim; i++) {
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);
    }

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 0; j < mesh->dim.xdim; j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 8 block */
void update_8(mesh_t *mesh, mesh_t *mesh_old, const cell_ops_t *cell_ops)
{
    int i, j;

    /* Update corner */
    cell_ops->cell_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, mesh->dim.xdim - 1, 2, 1);

    for (i = 0; i < (mesh->dim.xdim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);

    for (i = 0; i < (mesh->dim.ydim - 1); i++)
        cell_ops->cell_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_ops->cell_update_interior(mesh, mesh_old, i, j);
    }
}

void mesh_update(mesh_t *mesh, mesh_t *mesh_old, int block_type, const cell_ops_t *cell_ops)
{
    switch (block_type) {
        case 9:
            update_9(mesh, mesh_old, cell_ops);
            break;
        case 0:
            update_0(mesh, mesh_old, cell_ops);
            break;
        case 1:
            update_1(mesh, mesh_old, cell_ops);
            break;
        case 2:
            update_2(mesh, mesh_old, cell_ops);
            break;
        case 3:
            update_3(mesh, mesh_old, cell_ops);
            break;
        case 4:
            update_4(mesh, mesh_old, cell_ops);
            break;
        case 5:
            update_5(mesh, mesh_old, cell_ops);
            break;
        case 6:
            update_6(mesh, mesh_old, cell_ops);
            break;
        case 7:
            update_7(mesh, mesh_old, cell_ops);
            break;
        case 8:
            update_8(mesh, mesh_old, cell_ops);
            break;
    }
}

/* Checks for convergence at a specified cutoff. Returns 1 if relative error */
/* is less than the convergence cutoff, 0 otherwise */
int mesh_p_convergence_check(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff, int rank)
{
    int i, j;
    double num, denom, p_new, p_old, rel_error, global_num, global_denom;

    num = 0;
    denom = 0;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            p_new = mesh->cell[MESH_INDEX(i, j)].pressure;
            p_old = mesh_old->cell[MESH_INDEX(i, j)].pressure;
            num += pow(p_new - p_old, 2);
            denom += pow(p_new, 2);
        }
    }

    MPI_Reduce(&num, &global_num, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&denom, &global_denom, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        rel_error = sqrt(global_num / global_denom);
    }

    MPI_Bcast(&rel_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rel_error < conv_cutoff) {
        return 1;
    }
    else {
        return 0;
    }
}

/* Ensures the average pressure is 0 */
void mesh_p_impose_0_average(mesh_t *mesh, int rank)
{
    int N, i, j, k;
    double sum, global_sum, avg;

    N = mesh->dim.x_full_dim * mesh->dim.y_full_dim;
    sum = 0;
    global_sum = 0;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            sum += mesh->cell[MESH_INDEX(i, j)].pressure;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        avg = global_sum / N;
    }

    MPI_Bcast(&avg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            mesh->cell[MESH_INDEX(i, j)].pressure -= avg;

            for (k = 0; k < 4; k++) {
                mesh->cell[MESH_INDEX(i, j)].l_p[k] -= avg;
            }
        }
    }
}

/* Updates the robin conditions along the boundaries */
void mesh_update_robin(mesh_t *mesh, const cell_ops_t *cell_ops)
{
    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {
            cell_ops->cell_compute_robin(mesh, i, j);
        }
    }
}

/* Computes the velocity information on the mesh */
void mesh_compute_velocity(mesh_t *mesh)
{
    int i, j;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            cur_cell = &mesh->cell[MESH_INDEX(i, j)];

            cur_cell->velocity_y = (cur_cell->flux_p[0] - cur_cell->flux_p[2]) / 2;
            cur_cell->velocity_x = (cur_cell->flux_p[1] - cur_cell->flux_p[3]) / 2;
        }
    }
}

int mesh_pressure_iteration(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff,
    int block_type, int rank, send_vectors_t *send_vec, receive_vectors_t *rec_vec)
{
    int itr = 0;
    mesh_t *temp;

    mesh_compute_beta_A(mesh, &cell_p_ops);
    mesh_update_robin(mesh, &cell_p_ops);

    for (;;) {
        itr++;

        mesh_update(mesh, mesh_old, block_type, &cell_p_ops);

        if (mesh_p_convergence_check(mesh, mesh_old, conv_cutoff, rank)) {
            break;
        }

        mesh_p_impose_0_average(mesh, rank);
        mesh_update_robin(mesh, &cell_p_ops);
        mpi_comm(mesh, send_vec, rec_vec, block_type, rank);

        temp = mesh;
        mesh = mesh_old;
        mesh_old = temp;
    }

    return itr;
}

int mesh_diffusion_iteration(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff,
    int block_type, int rank, send_vectors_t *send_vec, receive_vectors_t *rec_vec)
{
    int itr = 0;
    mesh_t *temp;

    /* Computes diffusion D at all mesh points */
    for (int i = 0; i < mesh->dim.ydim; i++) {
        for (int j = 0; j < mesh->dim.xdim; j++) {

        }
    }

    return 0;
}

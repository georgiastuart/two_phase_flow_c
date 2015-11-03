#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "mpi.h"
#include "cell_functions.h"

/* Sets up the mesh with cell structs at each gridpoint */
mesh_t* mesh_init_mesh(dim_t dim, double *perm, double perm_scale, double perm_strength,
                        double *source, double c)
{
    int i, j;
    mesh_t *mesh;
    cell_t *cur_cell;

    /* Allocates memory for the mesh */
    mesh = malloc(sizeof(mesh_t));
    mesh->cell = malloc((dim.ydim + 2) * (dim.xdim + 2) * sizeof(cell_t));

    mesh->dim = dim;

    /* Sets cell permeability and sources */
    for (i = 0; i < (mesh->dim.ydim + 2); i++) {
        for (j = 0; j < (mesh->dim.xdim + 2); j++) {
            cur_cell = &mesh->cell[MESH_INDEX_INC_PAD(i, j)];
            cur_cell->perm = pow(10, -11) * exp(perm_strength * perm[MESH_INDEX_INC_PAD(i, j)]);
            cur_cell->source = source[MESH_INDEX_INC_PAD(i, j)];
        }
    }

    /* Computes beta and A at all mesh points */
    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            cell_p_compute_beta(mesh, i, j, c);
            cell_p_compute_A(mesh, i, j);
        }
    }

    return mesh;
}

/* One iteration of the algorithm over all mesh points */
/* 0 - up, 1 - right,  2 - down, 3 - left */
void iteration_9(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    /* Updates the corners */
    cell_p_update_corner(mesh, mesh_old, 0, 0, 0, 3);
    cell_p_update_corner(mesh, mesh_old, 0, mesh->dim.xdim - 1, 0, 1);
    cell_p_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, 0, 2, 3);
    cell_p_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, mesh->dim.xdim - 1, 2, 1);

    /* Updates the boundaries */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        cell_p_update_boundary(mesh, mesh_old, i, 0, 3);
        cell_p_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);
    }

    for (i = 1; i < (mesh->dim.xdim - 1); i++) {
        cell_p_update_boundary(mesh, mesh_old, 0, i, 0);
        cell_p_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);
    }

    /* Updates the interior cells */
    for (i = 1; i < (mesh->dim.ydim - 1); i++) {
        for (j = 1; j < (mesh->dim.xdim - 1); j++) {
            cell_p_update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* Iteration for a type 0 block */
void iteration_0(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    /* Update corner */
    cell_p_update_corner(mesh, mesh_old, 0, 0, 0, 3);

    for (i = 1; i < mesh->dim.xdim; i++)
        cell_p_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++)
        cell_p_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 1 block */
void iteration_1(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    for (i = 0; i < mesh->dim.xdim; i++)
        cell_p_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 2 block */
void iteration_2(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    /* Update corner */
    cell_p_update_corner(mesh, mesh_old, 0, mesh->dim.xdim - 1, 0, 1);

    for (i = 0; i < (mesh->dim.xdim - 1); i++)
        cell_p_update_boundary(mesh, mesh_old, 0, i, 0);

    for (i = 1; i < mesh->dim.ydim; i++)
        cell_p_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 1; i < mesh->dim.ydim; i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 3 block */
void iteration_3(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++)
        cell_p_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 4 block */
void iteration_4(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            cell_p_update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* Iteration for a type 5 block */
void iteration_5(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    for (i = 0; i < mesh->dim.ydim; i++)
        cell_p_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 6 block */
void iteration_6(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    /* Update corner */
    cell_p_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, 0, 2, 3);

    for (i = 1; i < mesh->dim.xdim; i++)
        cell_p_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);

    for (i = 0; i < (mesh->dim.ydim - 1); i++)
        cell_p_update_boundary(mesh, mesh_old, i, 0, 3);

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 1; j < mesh->dim.xdim; j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 7 block */
void iteration_7(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    for (i = 0; i < mesh->dim.xdim; i++) {
        cell_p_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);
    }

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 0; j < mesh->dim.xdim; j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

/* Iteration for a type 8 block */
void iteration_8(mesh_t *mesh, mesh_t *mesh_old)
{
    int i, j;

    /* Update corner */
    cell_p_update_corner(mesh, mesh_old, mesh->dim.ydim - 1, mesh->dim.xdim - 1, 2, 1);

    for (i = 0; i < (mesh->dim.xdim - 1); i++)
        cell_p_update_boundary(mesh, mesh_old, mesh->dim.ydim - 1, i, 2);

    for (i = 0; i < (mesh->dim.ydim - 1); i++)
        cell_p_update_boundary(mesh, mesh_old, i, mesh->dim.xdim - 1, 1);

    for (i = 0; i < (mesh->dim.ydim - 1); i++) {
        for (j = 0; j < (mesh->dim.xdim - 1); j++)
            cell_p_update_interior(mesh, mesh_old, i, j);
    }
}

void mesh_iteration(mesh_t *mesh, mesh_t *mesh_old, int block_type) {
    switch (block_type) {
        case 9:
            iteration_9(mesh, mesh_old);
            break;
        case 0:
            iteration_0(mesh, mesh_old);
            break;
        case 1:
            iteration_1(mesh, mesh_old);
            break;
        case 2:
            iteration_2(mesh, mesh_old);
            break;
        case 3:
            iteration_3(mesh, mesh_old);
            break;
        case 4:
            iteration_4(mesh, mesh_old);
            break;
        case 5:
            iteration_5(mesh, mesh_old);
            break;
        case 6:
            iteration_6(mesh, mesh_old);
            break;
        case 7:
            iteration_7(mesh, mesh_old);
            break;
        case 8:
            iteration_8(mesh, mesh_old);
            break;
    }
}

/* Checks for convergence at a specified cutoff. Returns 1 if relative error */
/* is less than the convergence cutoff, 0 otherwise */
int mesh_convergence_check(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff, int rank)
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
void mesh_impose_0_average(mesh_t *mesh, int rank)
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
                mesh->cell[MESH_INDEX(i, j)].l[k] -= avg;
            }
        }
    }
}

/* Updates the robin conditions along the boundaries */
void mesh_update_robin(mesh_t *mesh)
{
    int i, j, k;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.ydim; i++) {
        for (j = 0; j < mesh->dim.xdim; j++) {
            cur_cell = &mesh->cell[MESH_INDEX(i, j)];

            for (k = 0; k < 4; k++) {
                cur_cell->robin[k] = cur_cell->beta[k] * cur_cell->flux[k] + cur_cell->l[k];
            }
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

            cur_cell->velocity_y = (cur_cell->flux[0] - cur_cell->flux[2]) / 2;
            cur_cell->velocity_x = (cur_cell->flux[1] - cur_cell->flux[3]) / 2;
        }
    }
}

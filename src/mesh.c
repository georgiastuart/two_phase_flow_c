#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "cell_functions.h"

#define PERM_COEF pow(10, -11)

dim_t dim;

/* Sets up the mesh with cell structs at each gridpoint */
cell_t* init_mesh(double *perm, double perm_strength, double *source, double c)
{
    int i, j;
    cell_t *mesh, *cur_cell;

    /* Allocates memory for the mesh */
    mesh = malloc((dim.ydim + 2) * (dim.xdim + 2) * sizeof(cell_t));

    /* Sets cell permeability and sources */
    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            cur_cell = &mesh[MESH_INDEX(i, j)];
            cur_cell->perm = PERM_COEF * exp(perm_strength * perm[MESH_INDEX_NO_PAD(i, j)]);
            cur_cell->source = source[MESH_INDEX_NO_PAD(i, j)];
        }
    }

    /* Computes beta and A at all mesh points */
    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            compute_beta(mesh, i, j, c);
            compute_A(mesh, i, j);
        }
    }

    return mesh;
}

/* One iteration of the algorithm over all mesh points */
/* 0 - up, 1 - right,  2 - down, 3 - left */
void iteration(cell_t *mesh, cell_t *mesh_old)
{
    int i, j;

    /* Updates the corners */
    update_corner(mesh, mesh_old, 0, 0, 0, 3);
    update_corner(mesh, mesh_old, 0, dim.xdim - 1, 0, 1);
    update_corner(mesh, mesh_old, dim.ydim - 1, 0, 2, 3);
    update_corner(mesh, mesh_old, dim.ydim - 1, dim.xdim - 1, 2, 1);

    /* Updates the boundaries */
    for (i = 1; i < (dim.ydim - 1); i++) {
        update_boundary(mesh, mesh_old, i, 0, 3);
        update_boundary(mesh, mesh_old, i, dim.xdim - 1, 1);
    }

    for (i = 1; i < (dim.xdim - 1); i++) {
        update_boundary(mesh, mesh_old, 0, i, 0);
        update_boundary(mesh, mesh_old, dim.ydim - 1, i, 2);
    }

    /* Updates the interior cells */
    for (i = 1; i < (dim.ydim - 1); i++) {
        for (j = 1; j < (dim.xdim - 1); j++) {
            update_interior(mesh, mesh_old, i, j);
        }
    }
}

/* Checks for convergence at a specified cutoff. Returns 1 if relative error */
/* is less than the convergence cutoff, 0 otherwise */
int convergence_check(cell_t *mesh, cell_t *mesh_old, double conv_cutoff)
{
    int i, j;
    double num, denom, p_new, p_old, rel_error;

    num = 0;
    denom = 0;

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            p_new = mesh[MESH_INDEX(i, j)].pressure;
            p_old = mesh[MESH_INDEX(i, j)].pressure;
            num += pow(p_new - p_old, 2);
            denom += pow(p_new, 2);
        }
    }

    rel_error = sqrt(num / denom);

    if (rel_error < conv_cutoff) {
        return 1;
    }
    else {
        return 0;
    }
}

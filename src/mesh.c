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

/* Calculates an iteration for the pressure solution */

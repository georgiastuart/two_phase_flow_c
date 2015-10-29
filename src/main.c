#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"

#define DIM 8

int main(int argc, const char* argv[])
{
    double *perm, *source, perm_strength, beta_coef;
    int i, j;
    cell_t *mesh;

    dim.xdim = DIM;
    dim.ydim = DIM;
    perm_strength = 1;
    beta_coef = 1;

    /* Reads in the permeability and source fields */
    perm = read_file("perm_field_small.txt");
    source = read_file("src_field.txt");

    /* Initializes the mesh */
    mesh = init_mesh(perm, perm_strength, source, beta_coef);

    /* Frees memory used from the read-in permeability and source fields */
    free(perm);
    free(source);

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            printf("%e ", mesh[MESH_INDEX(i, j)].source);
        }
        printf("\n");
    }
    return 0;
}

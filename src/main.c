#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define DIM 8

int main(int argc, const char* argv[])
{
    double *perm, *source, perm_strength;
    dim_t dim;
    int i, j;
    cell_t *mesh;

    dim.xdim = DIM;
    dim.ydim = DIM;
    perm_strength = 1;

    /* Reads in the permeability and source fields */
    perm = read_file(dim, "perm_field_small.txt");
    source = read_file(dim, "src_field.txt");

    /* Initializes the mesh */
    mesh = setup_mesh(dim, perm, perm_strength, source);

    /* Frees memory used from the read-in permeability and source fields */
    free(perm);
    free(source);

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            printf("%e ", mesh[MESH_INDEX(i, j, dim.xdim)].source);
        }
        printf("\n");
    }
    return 0;
}

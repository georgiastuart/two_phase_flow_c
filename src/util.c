#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"

/* Sets up the mesh with cell structs at each gridpoint */
cell_t* setup_mesh(dim_t dim, double *perm, double *source)
{
    int i, j, xdim;
    cell_t *mesh, *cur_cell;
    mesh = malloc((dim.ydim + 2) * (dim.xdim + 2) * sizeof(cell_t));

    xdim = dim.xdim;

    for (i = 1; i < dim.ydim; i++) {
        for (j = 1; j < dim.xdim; j++) {
            cur_cell = &mesh[MESH_INDEX(i, j, xdim)];
            cur_cell->perm = perm[MESH_INDEX_NO_PAD(i, j, xdim)];
            cur_cell->source = source[MESH_INDEX_NO_PAD(i, j, xdim)];
        }
    }

    return mesh;
}

/* For reading permiability and source files */
double* read_file(dim_t dim, const char* file_name)
{
    FILE *fd;
    double *data;
    int i, j, xdim;

    data = malloc(dim.ydim * dim.xdim * sizeof(double));
    xdim = dim.xdim;

    if ((fd = fopen(file_name, "r")) == NULL) {
        printf("Cannot open file %s.\n", file_name);
        exit(1);
    }

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            fscanf(fd, "%lf", &data[MESH_INDEX_NO_PAD(i, j, xdim)]);
        }
    }

    fclose(fd);

    return data;
}

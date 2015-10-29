#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"

/* For reading permiability and source files */
double* read_file(const char* file_name)
{
    FILE *fd;
    double *data;
    int i, j;

    data = malloc(dim.ydim * dim.xdim * sizeof(double));

    if ((fd = fopen(file_name, "r")) == NULL) {
        printf("Cannot open file %s.\n", file_name);
        exit(1);
    }

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            fscanf(fd, "%lf", &data[MESH_INDEX_NO_PAD(i, j)]);
        }
    }

    fclose(fd);

    return data;
}

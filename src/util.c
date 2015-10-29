#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"

/* For reading permeability and source files */
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

/* Prints out the specified attribute to the terminal */
void print_attribute(cell_t *mesh, char *attribute)
{
    int i, j, k;
    if ( !strcmp(attribute, "beta") ) {
        for (k = 0; k < 4; k++) {
            switch (k) {
                case 0:
                    printf("Left Beta\n-------------------------------------\n");
                    break;
                case 1:
                    printf("Right Beta\n-------------------------------------\n");
                    break;
                case 2:
                    printf("Up Beta\n-------------------------------------\n");
                    break;
                case 3:
                    printf("Down Beta\n-------------------------------------\n");
                    break;
            }

            for (i = 0; i < dim.ydim; i++) {
                for (j = 0; j < dim.xdim; j++) {
                    printf("%e\t", mesh[MESH_INDEX(i, j)].beta[k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

    if ( !strcmp(attribute, "perm") ) {
        printf("PERMEABILITY\n------------------------------------------------\n");
        for (i = 0; i < dim.ydim; i++) {
            for (j = 0; j < dim.xdim; j++) {
                printf("%e\t", mesh[MESH_INDEX(i, j)].perm);
            }
            printf("\n");
        }
        printf("\n");
    }
}

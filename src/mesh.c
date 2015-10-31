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

/* Initializes the dimensions from the config file */
void init_dim(config_t *config)
{
    dim.xdim = config->xdim / config->num_subdomains_x;
    dim.ydim = config->ydim / config->num_subdomains_y;
    dim.xlen = config->xlen / config->num_subdomains_x;
    dim.ylen = config->ylen / config->num_subdomains_y;
    dim.h = dim.xlen / dim.xdim;
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
            p_old = mesh_old[MESH_INDEX(i, j)].pressure;
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

/* Ensures the average pressure is 0 */
void impose_0_average(cell_t *mesh)
{
    int N, i, j, k;
    double avg;

    N = dim.xdim * dim.ydim;
    avg = 0;

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            avg += mesh[MESH_INDEX(i, j)].pressure;
        }
    }

    avg /= N;

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            mesh[MESH_INDEX(i, j)].pressure -= avg;

            for (k = 0; k < 4; k++) {
                mesh[MESH_INDEX(i, j)].l[k] -= avg;
            }
        }
    }
}

/* Updates the robin conditions along the boundaries */
void update_robin(cell_t *mesh)
{
    int i, j, k;
    cell_t *cur_cell;

    for (i = 0; i < dim.ydim; i++) {
        for (j = 0; j < dim.xdim; j++) {
            cur_cell = &mesh[MESH_INDEX(i, j)];

            for (k = 0; k < 4; k++) {
                cur_cell->robin[k] = cur_cell->beta[k] * cur_cell->flux[k] + cur_cell->l[k];
            }
        }
    }
}

/* Prints out the specified attribute to the terminal */
void print_attribute(cell_t *mesh, char *attribute)
{
    int i, j, k;
    if ( !strcmp(attribute, "beta") ) {
        for (k = 0; k < 4; k++) {
            switch (k) {
                case 0:
                    printf("Up Beta\n-------------------------------------\n");
                    break;
                case 1:
                    printf("Right Beta\n-------------------------------------\n");
                    break;
                case 2:
                    printf("Down Beta\n-------------------------------------\n");
                    break;
                case 3:
                    printf("Left Beta\n-------------------------------------\n");
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

    if ( !strcmp(attribute, "pressure") ) {
        printf("PRESSURE\n----------------------------------------------------\n");
        for (i = 0; i < dim.ydim; i++) {
            for (j = 0; j < dim.xdim; j++) {
                printf("%e\t", mesh[MESH_INDEX(i, j)].pressure);
            }
            printf("\n");
        }
    }

    if ( !strcmp(attribute, "source") ) {
        printf("SOURCE\n-------------------------------------------------------\n");
        for (i = 0; i < dim.ydim; i++) {
            for (j = 0; j < dim.xdim; j++) {
                printf("%e\t", mesh[MESH_INDEX(i, j)].source);
            }
            printf("\n");
        }
    }

    if ( !strcmp(attribute, "flux") ) {
        for (k = 0; k < 4; k++) {
            switch (k) {
                case 0:
                    printf("Up Flux\n-------------------------------------\n");
                    break;
                case 1:
                    printf("Right Flux\n-------------------------------------\n");
                    break;
                case 2:
                    printf("Down Flux\n-------------------------------------\n");
                    break;
                case 3:
                    printf("Left Flux\n-------------------------------------\n");
                    break;
            }

            for (i = 0; i < dim.ydim; i++) {
                for (j = 0; j < dim.xdim; j++) {
                    printf("%e\t", mesh[MESH_INDEX(i, j)].flux[k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

    if ( !strcmp(attribute, "l") ) {
        for (k = 0; k < 4; k++) {
            switch (k) {
                case 0:
                    printf("Up l\n-------------------------------------\n");
                    break;
                case 1:
                    printf("Right l\n-------------------------------------\n");
                    break;
                case 2:
                    printf("Down l\n-------------------------------------\n");
                    break;
                case 3:
                    printf("Left l\n-------------------------------------\n");
                    break;
            }

            for (i = 0; i < dim.ydim; i++) {
                for (j = 0; j < dim.xdim; j++) {
                    printf("%e\t", mesh[MESH_INDEX(i, j)].l[k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

    if ( !strcmp(attribute, "robin") ) {
        for (k = 0; k < 4; k++) {
            switch (k) {
                case 0:
                    printf("Up Robin Conditions\n-------------------------------------\n");
                    break;
                case 1:
                    printf("Right Robin Conditions\n-------------------------------------\n");
                    break;
                case 2:
                    printf("Down Robin Conditions\n-------------------------------------\n");
                    break;
                case 3:
                    printf("Left Robin Conditions\n-------------------------------------\n");
                    break;
            }

            for (i = 0; i < dim.ydim; i++) {
                for (j = 0; j < dim.xdim; j++) {
                    printf("%e\t", mesh[MESH_INDEX(i, j)].robin[k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
}

/* Prints specified attribute to a file named "attribute".dat */
void print_attribute_to_file(cell_t *mesh, char *attribute)
{
    int i, j;
    FILE *fd;
    char name[100];

    sprintf(name, "output/%s.dat", attribute);

    if ((fd = fopen(name, "w")) == NULL) {
        fprintf(stderr, "Unable to open file %s\n", name);
        return;
    }

    if (!strcmp(attribute, "pressure")) {
        for (i = 0; i < dim.ydim; i++) {
            for (j = 0; j < dim.xdim; j++) {
                if (j == (dim.xdim - 1)) {
                    fprintf(fd, "%e\n", mesh[MESH_INDEX(i, j)].pressure);
                }
                else {
                    fprintf(fd, "%e,", mesh[MESH_INDEX(i, j)].pressure);
                }
            }
        }
    }

    fclose(fd);
}

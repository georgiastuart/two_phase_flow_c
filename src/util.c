#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

/* Necessary function for the config reader */
static int config_helper(void *config, const char *section, const char *name,
                            const char *value)
{
    config_t *pconfig = (config_t*)config;
    double d;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0

    if (MATCH("dimensions", "xdim")) {
        pconfig->xdim = atoi(value);
    } else if (MATCH("dimensions", "ydim")) {
        pconfig->ydim = atoi(value);
    } else if (MATCH("dimensions", "xlen")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->xlen = d;
    } else if (MATCH("dimensions", "ylen")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->ylen = d;
    } else if (MATCH("files","perm_file")) {
        strcpy(pconfig->perm_file, strdup(value));
    } else if (MATCH("files","src_file")) {
        strcpy(pconfig->src_file, strdup(value));
    } else if (MATCH("other", "perm_scale")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->perm_scale = d;
    } else if (MATCH("other", "perm_strength")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->perm_strength = d;
    } else if (MATCH("other", "conv_cutoff")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->conv_cutoff = d;
    } else if (MATCH("other", "beta_coef")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->beta_coef = d;
    } else {
        /* Unknown Section or Name */
        return 0;
    }

    return 1;
}

/* Reads the configuration file */
int read_config(const char* file_name, config_t *config)
{

    if (ini_parse(file_name, config_helper, config) < 0) {
        printf("Can't load %s\n", file_name);
        return 0;
    }
    printf("Config loaded from %s\n", file_name);

    return 1;
}

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

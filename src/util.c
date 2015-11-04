#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"

#define INDEX(y, x, xdim) (y * xdim + x)
#define INDEX_PAD(y, x, xdim) ((y + 1) * (xdim + 2) + (x + 1))

/* Initializes the dimensions from the config file */
void init_dim(config_t *config, dim_t *dim)
{
    dim->xdim = config->xdim / config->num_subdomains_x;
    dim->ydim = config->ydim / config->num_subdomains_y;
    dim->x_full_dim = config->xdim;
    dim->y_full_dim = config->ydim;
    dim->xlen = config->xlen / config->num_subdomains_x;
    dim->ylen = config->ylen / config->num_subdomains_y;
    dim->num_subdomains_x = config->num_subdomains_x;
    dim->num_subdomains_y = config->num_subdomains_y;
    dim->h = dim->xlen / dim->xdim;
}

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
        strcpy(pconfig->perm_file, value);
    } else if (MATCH("files","src_file")) {
        strcpy(pconfig->src_file, value);
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
    } else if (MATCH("other", "porosity")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->porosity = d;
    } else if (MATCH("other", "sat_rel_o")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->sat_rel_o = d;
    } else if (MATCH("other", "sat_rel_w")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->sat_rel_w = d;
    } else if (MATCH("other", "visc_o")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->visc_o = d;
    } else if (MATCH("other", "visc_w")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->visc_w = d;
    } else if (MATCH("other", "eta")) {
        sscanf(strdup(value), "%lf", &d);
        pconfig->eta = d;
    } else if (MATCH("mpi", "num_processes")) {
        pconfig->num_processes = atoi(value);
    } else if (MATCH("mpi", "num_subdomains_x")) {
        pconfig->num_subdomains_x = atoi(value);
    } else if (MATCH("mpi", "num_subdomains_y")) {
        pconfig->num_subdomains_y = atoi(value);
    } else if (MATCH("out_files", "pressure_out")) {
        strcpy(pconfig->pressure_out, value);
    } else if (MATCH("out_files", "velocity_y_out")) {
        strcpy(pconfig->velocity_y_out, value);
    } else if (MATCH("out_files", "velocity_x_out")) {
        strcpy(pconfig->velocity_x_out, value);
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

    if ((config->num_subdomains_x * config->num_subdomains_y) != config->num_processes) {
        printf("Number of subdomains does not match number of processes\n");
        return 0;
    }

    printf("Config loaded from %s\n", file_name);

    return 1;
}

/* For reading permeability and source files */
double* read_file(const char* file_name, int ydim, int xdim)
{
    FILE *fd;
    double *data;
    int i, j;

    data = malloc((ydim + 2) * (xdim + 2) * sizeof(double));

    if ((fd = fopen(file_name, "r")) == NULL) {
        printf("Cannot open file %s.\n", file_name);
        exit(1);
    }

    for (i = 0; i < (ydim + 2); i++) {
        for (j = 0; j < (xdim + 2); j++) {
            fscanf(fd, "%lf", &data[INDEX(i, j, (xdim + 2))]);
        }
    }

    fclose(fd);

    return data;
}

/* For reading and padding the original permeability and source files */
double* read_file_pad(const char* file_name, int ydim, int xdim)
{
    FILE *fd;
    double *data;
    int i, j;

    data = malloc((ydim + 2) * (xdim + 2) * sizeof(double));

    if ((fd = fopen(file_name, "r")) == NULL) {
        printf("Cannot open file %s.\n", file_name);
        exit(1);
    }

    for (i = 0; i < ydim; i++) {
        for (j = 0; j < xdim; j++) {
            fscanf(fd, "%lf", &data[INDEX_PAD(i, j, xdim)]);
        }
    }

    fclose(fd);

    return data;
}

void setup_files(const char* file_name, int ydim, int xdim, int num_subdomains_y,
                    int num_subdomains_x, int size, const char* mode)
{
    double *temp;
    FILE *read, *write;
    char name[100];
    int i, j, k, xdim_per_block, ydim_per_block, x_block_loc, y_block_loc;
    double *val;

    xdim_per_block = xdim / num_subdomains_x;
    ydim_per_block = ydim / num_subdomains_y;

    if ((read = fopen(file_name, "r")) == NULL) {
        printf("Cannot open file %s.\n", file_name);
        exit(1);
    }

    temp = read_file_pad(file_name, ydim, xdim);

    for (k = 0; k < size; k++) {
        x_block_loc = k % num_subdomains_x;
        y_block_loc = (int)(k / num_subdomains_y);

        sprintf(name, "input/%s.%d", mode, k);

        if ((write = fopen(name, "w")) == NULL) {
            printf("Cannot open file %s.\n", name);
            exit(1);
        }

        for (i = 0; i < (ydim_per_block + 2); i++) {
            for (j = 0; j < (xdim_per_block + 2); j++) {
                val = &temp[INDEX((i + y_block_loc * ydim_per_block),
                                    (j + x_block_loc * xdim_per_block), (xdim + 2))];
                fprintf(write, "%e\n", *val);
            }
        }
        fclose(write);
    }

    fclose(read);
    free(temp);
}

/* Prints out the specified attribute to the terminal */
void print_attribute(mesh_t *mesh, char *attribute)
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

            for (i = 0; i < mesh->dim.ydim; i++) {
                for (j = 0; j < mesh->dim.xdim; j++) {
                    printf("%e\t", mesh->cell[MESH_INDEX(i, j)].beta[k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }

    if ( !strcmp(attribute, "perm") ) {
        printf("PERMEABILITY\n------------------------------------------------\n");
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                printf("%e\t", mesh->cell[MESH_INDEX(i, j)].perm);
            }
            printf("\n");
        }
        printf("\n");
    }

    if ( !strcmp(attribute, "pressure") ) {
        printf("PRESSURE\n----------------------------------------------------\n");
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                printf("%e\t", mesh->cell[MESH_INDEX(i, j)].pressure);
            }
            printf("\n");
        }
    }

    if ( !strcmp(attribute, "source") ) {
        printf("SOURCE\n-------------------------------------------------------\n");
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                printf("%e\t", mesh->cell[MESH_INDEX(i, j)].source);
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

            for (i = 0; i < mesh->dim.ydim; i++) {
                for (j = 0; j < mesh->dim.xdim; j++) {
                    printf("%e\t", mesh->cell[MESH_INDEX(i, j)].flux_p[k]);
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

            for (i = 0; i < mesh->dim.ydim; i++) {
                for (j = 0; j < mesh->dim.xdim; j++) {
                    printf("%e\t", mesh->cell[MESH_INDEX(i, j)].l_p[k]);
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

            for (i = 0; i < mesh->dim.ydim; i++) {
                for (j = 0; j < mesh->dim.xdim; j++) {
                    printf("%e\t", mesh->cell[MESH_INDEX(i, j)].robin[k]);
                }
                printf("\n");
            }
            printf("\n");
        }
    }
}

/* Prints specified attribute to a file named "attribute".dat */
void print_attribute_to_file(mesh_t *mesh, char *attribute)
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
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                if (j == (mesh->dim.xdim - 1)) {
                    fprintf(fd, "%e\n", mesh->cell[MESH_INDEX(i, j)].pressure);
                }
                else {
                    fprintf(fd, "%e,", mesh->cell[MESH_INDEX(i, j)].pressure);
                }
            }
        }
    }

    fclose(fd);
}

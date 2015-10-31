#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"

#define INDEX(y, x, xdim) (y * xdim + x)

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
    } else if (MATCH("mpi", "num_processes")) {
        pconfig->num_processes = atoi(value);
    } else if (MATCH("mpi", "num_subdomains_x")) {
        pconfig->num_subdomains_x = atoi(value);
    } else if (MATCH("mpi", "num_subdomains_y")) {
        pconfig->num_subdomains_y = atoi(value);
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

    data = malloc(ydim * xdim * sizeof(double));

    if ((fd = fopen(file_name, "r")) == NULL) {
        printf("Cannot open file %s.\n", file_name);
        exit(1);
    }

    for (i = 0; i < ydim; i++) {
        for (j = 0; j < xdim; j++) {
            fscanf(fd, "%lf", &data[i * xdim + j]);
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

    temp = read_file(file_name, ydim, xdim);

    for (k = 0; k < size; k++) {
        x_block_loc = k % num_subdomains_x;
        y_block_loc = (int)(k / num_subdomains_y);

        sprintf(name, "input/%s.%d", mode, k);

        if ((write = fopen(name, "w")) == NULL) {
            printf("Cannot open file %s.\n", name);
            exit(1);
        }

        for (i = 0; i < ydim_per_block; i++) {
            for (j = 0; j < xdim_per_block; j++) {
                val = &temp[INDEX((i + y_block_loc * ydim_per_block),
                                    (j + x_block_loc * xdim_per_block), xdim)];
                fprintf(write, "%e\n", *val);
            }
        }
        fclose(write);
    }

    fclose(read);
}

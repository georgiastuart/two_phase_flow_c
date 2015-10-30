#ifndef H_UTIL
#define H_UTIL

#include "mesh.h"
#include "ini.h"

typedef struct {
    int xdim, ydim;
    double xlen, ylen;
    char perm_file[100], src_file[100];
    double perm_scale, perm_strength, conv_cutoff, beta_coef;
} config_t;

int read_config(const char* file_name, config_t *config);
double* read_file(const char* file_name);
void print_attribute(cell_t *mesh, char *attribute);
void print_attribute_to_file(cell_t *mesh, char *attribute);

#endif /* H_UTIL */

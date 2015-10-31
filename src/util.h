#ifndef H_UTIL
#define H_UTIL

#include "ini.h"

typedef struct config
{
    int xdim, ydim;
    double xlen, ylen;
    char perm_file[100], src_file[100];
    double perm_scale, perm_strength, conv_cutoff, beta_coef;
    int num_processes, num_subdomains_x, num_subdomains_y;
} config_t;

int read_config(const char* file_name, config_t *config);
double* read_file(const char* file_name, int ydim, int xdim);
void setup_files(const char* file_name, int ydim, int xdim, int num_subdomains_y,
                    int num_subdomains_x, int size, const char* mode);

#endif /* H_UTIL */

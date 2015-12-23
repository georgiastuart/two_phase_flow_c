#ifndef H_UTIL
#define H_UTIL

#include "inih/ini.h"
#include "mpi.h"

struct mesh;
typedef struct mesh mesh_t;

typedef struct dim
{
    int xdim, ydim, num_subdomains_x, num_subdomains_y, x_full_dim, y_full_dim;
    int num_ts;
    double h, xlen, ylen, dt, dt_transport;
} dim_t;

typedef struct config
{
    int xdim, ydim, time_steps;
    double xlen, ylen, dt;
    char perm_file[100], src_file[100], pressure_out[100], velocity_y_out[100];
    char velocity_x_out[100], saturation_out[100], sat_file[100];
    double perm_scale, perm_strength, conv_cutoff, beta_coef;
    int linearity;
    int num_processes, num_subdomains_x, num_subdomains_y;
    double porosity, visc_o, visc_w, sat_rel_o, sat_rel_w, eta;
} config_t;

void init_dim(config_t *config, dim_t *dim);
int read_config(const char* file_name, config_t *config);
double* read_file(const char* file_name, int ydim, int xdim);
double* read_file_pad(const char* file_name, int ydim, int xdim);
void setup_files(const char* file_name, int ydim, int xdim, int num_subdomains_y,
                    int num_subdomains_x, int size, const char* mode);
void print_attribute(mesh_t *mesh, char *attribute);
void progress_bar(int cur_time, int total_time);

#endif /* H_UTIL */

#ifndef H_MESH
#define H_MESH

#define MESH_INDEX(y, x) ((y + 1) * (dim.xdim + 2) + (x + 1))
#define MESH_INDEX_NO_PAD(y, x) (y * dim.xdim + x)
#define MESH_INDEX_INC_PAD(y, x) (y * (dim.xdim + 2) + x)

#include "util.h"
#include "mpi.h"

typedef struct dim
{
    int xdim, ydim, num_subdomains_x, num_subdomains_y, x_full_dim, y_full_dim;
    double h, xlen, ylen;
} dim_t;

typedef struct cell
{
    /* Values that live at the center of each cell */
    double perm, pressure, source;

    /* Values that live along the edges */
    double flux[4], l[4], beta[4], robin[4], A[4];
} cell_t;

extern dim_t dim;

cell_t* init_mesh(double *perm, double perm_strength, double *source, double c);
void init_dim(config_t *config);
void iteration(cell_t *mesh, cell_t *mesh_old, int block_type);
int convergence_check(cell_t *mesh, cell_t *mesh_old, double conv_cutoff, int rank);
void impose_0_average(cell_t *mesh, int rank);
void update_robin(cell_t *mesh);
void print_attribute(cell_t *mesh, char *attribute);
void print_attribute_to_file(cell_t *mesh, char *attribute);

#endif /* H_MESH */

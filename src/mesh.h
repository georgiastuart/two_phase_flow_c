#ifndef H_MESH
#define H_MESH

#define MESH_INDEX(y, x) ((y + 1) * (mesh->dim.xdim + 2) + (x + 1))
#define MESH_INDEX_NO_PAD(y, x) (y * mesh->dim.xdim + x)
#define MESH_INDEX_INC_PAD(y, x) (y * (mesh->dim.xdim + 2) + x)

#include "util.h"

typedef struct cell
{
    /* Values that live at the center of each cell */
    double perm, pressure, source;

    /* Values that live along the edges */
    double flux[4], l[4], beta[4], robin[4], A[4];
} cell_t;

typedef struct mesh
{
    cell_t *cell;
    dim_t dim;
} mesh_t;

mesh_t* init_mesh(dim_t dim, double *perm, double perm_strength, double *source, double c);
void iteration(mesh_t *mesh, mesh_t *mesh_old, int block_type);
int convergence_check(mesh_t *mesh, mesh_t *mesh_old, double conv_cutoff, int rank);
void impose_0_average(mesh_t *mesh, int rank);
void update_robin(mesh_t *mesh);
void print_attribute(mesh_t *mesh, char *attribute);
void print_attribute_to_file(mesh_t *mesh, char *attribute);

#endif /* H_MESH */

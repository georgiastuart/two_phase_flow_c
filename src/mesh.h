#ifndef H_MESH
#define H_MESH

#define MESH_INDEX(y, x) ((y + 1) * (dim.xdim + 2) + (x + 1))
#define MESH_INDEX_NO_PAD(y, x) (y * dim.xdim + x)

typedef struct dim
{
    int xdim, ydim;
    double h, xphysdim, yphysdim;
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
void iteration(cell_t *mesh, cell_t *mesh_old);
int convergence_check(cell_t *mesh, cell_t *mesh_old, double conv_cutoff);
void impose_0_average(cell_t *mesh);
void update_robin(cell_t *mesh);

#endif /* H_MESH */

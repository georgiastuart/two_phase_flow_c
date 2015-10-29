#ifndef H_MESH
#define H_MESH

#define MESH_INDEX(y, x) ((y + 1) * (dim.xdim + 2) + (x + 1))
#define MESH_INDEX_NO_PAD(y, x) (y * dim.xdim + x)

typedef struct dim
{
    int xdim, ydim;
    double h, x_left_bound, x_right_bound, y_up_bound, y_low_bound;
} dim_t;

typedef struct cell
{
    double perm, pressure, source;
    double q[4], l[4], beta[4];
} cell_t;

extern dim_t dim;

cell_t* init_mesh(double *perm, double perm_strength, double *source, double c);

#endif /* H_MESH */

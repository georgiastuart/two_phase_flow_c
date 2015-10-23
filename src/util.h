#ifndef H_UTIL
#define H_UTIL

#define MESH_INDEX(y, x, xdim) (y * (xdim + 2) + x)
#define MESH_INDEX_NO_PAD(y, x, xdim) (y * xdim + x)

typedef struct dim
{
    int xdim, ydim;
} dim_t;

typedef struct edge
{
    double left, right, up, down;
} edge_t;

typedef struct cell
{
    double perm, pressure, source;
    edge_t q, l, beta;
} cell_t;

cell_t* setup_mesh(dim_t dim, double *perm, double *source);
double* read_file(dim_t dim, const char* file_name);

#endif /* H_UTIL */

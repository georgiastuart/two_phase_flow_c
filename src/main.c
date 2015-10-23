#include <stdio.h>
#include "util.h"

#define DIM 64

int main(int argc, const char* argv[])
{
    double *perm, *source;
    FILE *fd;
    dim_t dim;
    int i, j;
    cell_t mesh*;

    dim.xdim = DIM;
    dim.ydim = DIM;

    perm = read_file(dim, "perm_field.txt");
    source = read_file(dim, "src_field.txt");

    mesh = setup_mesh(dim, perm, source);

    free(prem);
    free(source);
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"

#define DIM 8
#define PHYSDIM 25600

int main(int argc, const char* argv[])
{
    double *perm, *source, perm_strength, beta_coef;
    cell_t *mesh, *mesh_old;

    dim.xdim = DIM;
    dim.ydim = DIM;
    dim.xphysdim = PHYSDIM;
    dim.yphysdim = PHYSDIM;
    dim.h = PHYSDIM / DIM;
    perm_strength = 1;
    beta_coef = 1;

    /* Reads in the permeability and source fields */
    perm = read_file("perm_field_small.txt");
    source = read_file("src_field_small.txt");

    /* Initializes the meshes */
    mesh = init_mesh(perm, perm_strength, source, beta_coef);
    mesh_old = init_mesh(perm, perm_strength, source, beta_coef);

    /* Frees memory used from the read-in permeability and source fields */
    free(perm);
    free(source);

    print_attribute(mesh, "pressure");

    iteration(mesh, mesh_old);

    print_attribute(mesh, "pressure");

    print_attribute(mesh, "source");

    return 0;
}

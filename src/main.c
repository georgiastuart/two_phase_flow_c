#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"

int main(int argc, const char* argv[])
{
    double *perm, *source;
    cell_t *mesh, *mesh_old, *temp;
    config_t config;

    if (!read_config("input/config.ini", &config))
        return 1;

    dim.xdim = config.xdim;
    dim.ydim = config.ydim;
    dim.xlen = config.xlen;
    dim.ylen = config.ylen;
    dim.h = dim.xlen / dim.xdim;

    /* Reads in the permeability and source fields */
    perm = read_file(config.perm_file);
    source = read_file(config.src_file);

    /* Initializes the meshes */
    mesh = init_mesh(perm, config.perm_strength, source, config.beta_coef);
    mesh_old = init_mesh(perm, config.perm_strength, source, config.beta_coef);

    /* Frees memory used from the read-in permeability and source fields */
    free(perm);
    free(source);

    int itr;
    itr = 0;
    for (;;) {
        iteration(mesh, mesh_old);
        if ( convergence_check(mesh, mesh_old, config.conv_cutoff) ) {
            break;
        }

        impose_0_average(mesh);
        update_robin(mesh);

        temp = mesh;
        mesh = mesh_old;
        mesh_old = temp;

        itr++;
    }

    printf("Finished after %d iterations.\n", itr + 1);
    print_attribute(mesh, "pressure");
    print_attribute_to_file(mesh, "pressure");

    return 0;
}

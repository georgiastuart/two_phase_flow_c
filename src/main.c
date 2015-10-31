#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"
#include "mpi_util.h"
#include "mpi.h"

int main(int argc, char* argv[])
{
    double *perm, *source;
    cell_t *mesh, *mesh_old, *temp;
    config_t config;
    int rank, size;
    MPI_Datatype mpi_config_t;

    /* Initializes MPI and creates the config datatype */
    mpi_setup(&argc, &argv, &rank, &size, &mpi_config_t);

    /* Reads in config file */
    if (rank == 0) {
        if (!read_config("input/config.ini", &config))
            return 1;
    }

    MPI_Bcast(&config, 1, mpi_config_t, 0, MPI_COMM_WORLD);

    printf("%d\n", config.xdim);

    init_dim(&config);

    /* Reads in the permeability and source fields */
    perm = read_file(config.perm_file, dim.ydim, dim.xdim);
    source = read_file(config.src_file, dim.ydim, dim.xdim);

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

    if (rank == 0) {
        printf("Finished after %d iterations.\n", itr + 1);
        print_attribute(mesh, "pressure");
        print_attribute_to_file(mesh, "pressure");
    }

    mpi_shutdown(&mpi_config_t);

    return 0;
}

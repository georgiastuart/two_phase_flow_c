#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"
#include "mpi_util.h"

int main(int argc, char* argv[])
{
    double *perm, *source;
    mesh_t *mesh, *mesh_old, *temp;
    dim_t dim;
    config_t config;
    int rank, size, block_type;
    MPI_Datatype mpi_config_t;
    receive_vectors_t rec_vec;
    send_vectors_t send_vec;
    char perm_file[100], src_file[100];

    /* Initializes MPI and creates the config datatype */
    mpi_setup(&argc, &argv, &rank, &size, &mpi_config_t);

    /* Reads in config file */
    if (rank == 0) {
        if (!read_config("input/config.ini", &config))
            return 1;
    }

    MPI_Bcast(&config, 1, mpi_config_t, 0, MPI_COMM_WORLD);

    /* Sets up files for perm and source */
    if (rank == 0) {
        setup_files(config.perm_file, config.ydim, config.xdim, config.num_subdomains_y,
                    config.num_subdomains_x, size, "perm");
        setup_files(config.src_file, config.ydim, config.xdim, config.num_subdomains_y,
                    config.num_subdomains_x, size, "src");
    }

    /* Gets the type of subdomain */
    block_type = get_block_type(rank, config.num_subdomains_y, config.num_subdomains_x);

    printf("%d, %d, %d\n", config.xdim, block_type, rank);

    /* Changes perm and source file name to correct split file */
    sprintf(perm_file, "input/perm.%d", rank);
    sprintf(src_file, "input/src.%d", rank);

    /* Initializes the global dimension struct */
    init_dim(&config, &dim);

    /* Reads in the permeability and source fields */
    perm = read_file(perm_file, dim.ydim, dim.xdim);
    source = read_file(src_file, dim.ydim, dim.xdim);

    /* Initializes the meshes */
    mesh = init_mesh(dim, perm, config.perm_strength, source, config.beta_coef);
    mesh_old = init_mesh(dim, perm, config.perm_strength, source, config.beta_coef);

    /* Initializes send and receive structs for MPI */
    init_send_receive(mesh, &send_vec, &rec_vec);

    /* Frees memory used from the read-in permeability and source fields */
    free(perm);
    free(source);

    int itr;
    itr = 0;
    for (;;) {
        iteration(mesh, mesh_old, block_type);
        if ( convergence_check(mesh, mesh_old, config.conv_cutoff, rank) ) {
            break;
        }

        impose_0_average(mesh, rank);

        update_robin(mesh);

        comm(mesh, &send_vec, &rec_vec, block_type, rank);

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

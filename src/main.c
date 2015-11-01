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
    int rank, size, block_type, is_master;
    MPI_Datatype mpi_config_t;
    receive_vectors_t rec_vec;
    send_vectors_t send_vec;
    char perm_file[100], src_file[100];
    double t1, t2;
    int i, j;

    /* Initializes MPI and creates the config datatype */
    mpi_setup(&argc, &argv, &rank, &size, &mpi_config_t);

    /* Identifies the master process */
    is_master = (rank == 0);

    if (is_master) {
        printf("Setting up MPI and mesh...\n");
        t1 = MPI_Wtime();

        /* Reads in config file */
        if (!read_config("input/config.ini", &config))
            return 1;
    }

    /* Broadcasts config file */
    MPI_Bcast(&config, 1, mpi_config_t, 0, MPI_COMM_WORLD);

    /* Sets up the permeability and source parameters and sends to processes */
    mpi_setup_parameters(&config, 0, size, is_master, &perm);
    mpi_setup_parameters(&config, 1, size, is_master, &source);

    /* Gets the type of subdomain */
    block_type = mpi_get_block_type(rank, config.num_subdomains_y, config.num_subdomains_x);

    /* Initializes the global dimension struct */
    init_dim(&config, &dim);

    /* Initializes the meshes */
    mesh = mesh_init_mesh(dim, perm, config.perm_scale, config.perm_strength,
                            source, config.beta_coef);
    mesh_old = mesh_init_mesh(dim, perm, config.perm_scale, config.perm_strength,
                            source, config.beta_coef);

    /* Initializes send and receive structs for MPI */
    mpi_init_send_receive(mesh, &send_vec, &rec_vec);

    /* Frees memory used from the read-in permeability and source fields */
    free(perm);
    free(source);

    if (is_master) {
        t2 = MPI_Wtime();
        printf("Setup finished after %f seconds\n", t2 - t1);
        printf("Running iterations...\n");
    }

    int itr;
    itr = 0;
    for (;;) {
        mesh_iteration(mesh, mesh_old, block_type);
        if ( mesh_convergence_check(mesh, mesh_old, config.conv_cutoff, rank) ) {
            break;
        }

        mesh_impose_0_average(mesh, rank);

        mesh_update_robin(mesh);

        mpi_comm(mesh, &send_vec, &rec_vec, block_type, rank);

        temp = mesh;
        mesh = mesh_old;
        mesh_old = temp;

        itr++;
    }

    if (rank == 0) {
        t1 = MPI_Wtime();
        printf("Finished after %f seconds and %d iterations.\n", t1 - t2, itr + 1);
        /*print_attribute(mesh, "pressure");*/
        print_attribute_to_file(mesh, "pressure");
    }

    mpi_shutdown(&mpi_config_t);

    return 0;
}

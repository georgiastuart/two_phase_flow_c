#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "mesh.h"
#include "cell_functions.h"
#include "mpi_util.h"

#define NDIMS 2

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
    double t1, t2;

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
    printf("Before setup parameters\n");
    mpi_setup_parameters(&config, 0, size, is_master, &perm);
    printf("Setup perm\n");
    mpi_setup_parameters(&config, 1, size, is_master, &source);
    printf("Setup parameters\n");

    /* Gets the type of subdomain */
    block_type = mpi_get_block_type(rank, config.num_subdomains_y, config.num_subdomains_x);

    /* Initializes the global dimension struct */
    init_dim(&config, &dim);

    /* Initializes the meshes */
    mesh = mesh_init_mesh(dim, perm, source, &config);
    mesh_old = mesh_init_mesh(dim, perm, source, &config);

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

    /* Sets up mesh for diffusion test */
    setup_diffusion_test(mesh);
    setup_diffusion_test(mesh_old);

    int itr;

    for (int i = 0; i < config.time_steps; i++) {
        itr = mesh_diffusion_iteration(mesh, mesh_old, config.conv_cutoff, block_type,
                                       rank, &send_vec, &rec_vec);

        update_saturation_time(mesh, mesh_old);

        temp = mesh;
        mesh = mesh_old;
        mesh_old = temp;
    }

    // /* Iteration of the pressure problem */
    // int itr = mesh_pressure_iteration(mesh, mesh_old, config.conv_cutoff,
    //             block_type, rank, &send_vec, &rec_vec);
    //
    //
    if (is_master) {
        t1 = MPI_Wtime();
        printf("Finished after %f seconds and %d iterations.\n", t1 - t2, itr);
    }

    // /* Computes velocity */
    // mesh_compute_velocity(mesh);

    // /* Writes out data to binaries */
    // write_data(mesh, &config, size, rank, "pressure");
    // write_data(mesh, &config, size, rank, "velocity_y");
    // write_data(mesh, &config, size, rank, "velocity_x");
    write_data(mesh, &config, size, rank, "saturation");


    /* Shuts down MPI and frees memory */
    mpi_shutdown(&mpi_config_t);
    free(mesh);
    free(mesh_old);

    return 0;
}

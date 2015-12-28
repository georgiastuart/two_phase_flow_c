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
    double *perm, *source, *sat;
    mesh_t *mesh, *mesh_old;
    dim_t dim;
    config_t config;
    int rank, size, block_type, is_master;
    MPI_Datatype mpi_config_t;
    receive_vectors_t rec_vec;
    send_vectors_t send_vec;
    double t1, t2;
    production_wells_t prod_wells;

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
    mpi_setup_parameters(&config, 2, size, is_master, &sat);

    /* Gets the type of subdomain */
    block_type = mpi_get_block_type(rank, config.num_subdomains_y, config.num_subdomains_x);

    /* Initializes the global dimension struct */
    init_dim(&config, &dim);

    /* Sets linear or nonlinear parameters */
    set_linearity(config.linearity);

    /* Initializes the meshes */
    mesh = mesh_init_mesh(dim, perm, source, sat, &config);
    mesh_old = mesh_init_mesh(dim, perm, source, sat, &config);

    /* Initializes send and receive structs for MPI */
    mpi_init_send_receive(mesh, &send_vec, &rec_vec);

    /* Frees memory used from the read-in permeability and source fields */
    free(perm);
    free(source);

    /* Initializes Production Wells */
    prod_wells = init_production_wells(mesh);

    // printf("Num prod: %d\n", prod_wells.num_wells);

    if (is_master) {
        t2 = MPI_Wtime();
        printf("Setup finished after %f seconds\n", t2 - t1);
        printf("Running iterations...\n");
    }

    /* Two phase flow time steps */
    int itr;
    for (int i = 1; i < (config.time_steps + 1); i++) {
        if (is_master)
            progress_bar(i, config.time_steps);

        itr = mesh_pressure_iteration(mesh, mesh_old, config.conv_cutoff,
                    block_type, rank, &send_vec, &rec_vec);

        // if (is_master) {
        //     t1 = MPI_Wtime();
        //     printf("Time %d: Pressure finished after %f seconds and %d iterations.\n",
        //                 i, t1 - t2, itr);
        // }

        itr = mesh_transport_iteration(mesh, mesh_old, block_type, rank, &send_vec, &rec_vec);

        // if (is_master) {
        //     t2 = MPI_Wtime();
        //     printf("Time %d: Transport finished after %f seconds and %d iterations.\n",
        //                 i, t2 - t1, itr);
        // }

        if (!config.linearity) {
            itr = mesh_diffusion_iteration(mesh, mesh_old, config.conv_cutoff, block_type,
                                               rank, &send_vec, &rec_vec);
        }


        // if (is_master) {
        //     t1 = MPI_Wtime();
        //     printf("Time %d: Diffusion finished after %f seconds and %d iterations.\n",
        //                 i, t1 - t2, itr);
        // }

        record_production_wells(&prod_wells, mesh, i);
    }

    if (is_master) {
        t1 = MPI_Wtime();
        double time = t1 - t2;
        int mins = (int) (time / 60.0);
        double secs = time - 60.0 * mins;
        printf("Finished after %d minutes and %.02f seconds\n", mins, secs);
    }

    /* Writes out data to binaries */
    write_data(mesh, &config, size, rank, "pressure");
    write_data(mesh, &config, size, rank, "velocity_y");
    write_data(mesh, &config, size, rank, "velocity_x");
    write_data(mesh, &config, size, rank, "saturation");
    write_prod_well_data(&prod_wells, &config, size, rank);


    /* Shuts down MPI and frees memory */
    mpi_shutdown(&mpi_config_t);
    free(mesh);
    free(mesh_old);

    return 0;
}

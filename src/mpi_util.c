#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi_util.h"

void mpi_setup(int *argc, char ***argv, int *rank, int *size, MPI_Datatype *mpi_config_t)
{
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, size);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);

    /* Create type for config struct */
    const int nitems = 13;
    int blocklengths[13] = {1, 1, 1, 1, 100, 100, 1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype types[13] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE,
                                MPI_CHAR, MPI_CHAR,
                                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                MPI_INT, MPI_INT, MPI_INT};
    MPI_Aint offsets[13];

    offsets[0] = offsetof(config_t, xdim);
    offsets[1] = offsetof(config_t, ydim);
    offsets[2] = offsetof(config_t, xlen);
    offsets[3] = offsetof(config_t, ylen);
    offsets[4] = offsetof(config_t, perm_file);
    offsets[5] = offsetof(config_t, src_file);
    offsets[6] = offsetof(config_t, perm_scale);
    offsets[7] = offsetof(config_t, perm_strength);
    offsets[8] = offsetof(config_t, conv_cutoff);
    offsets[9] = offsetof(config_t, beta_coef);
    offsets[10] = offsetof(config_t, num_processes);
    offsets[11] = offsetof(config_t, num_subdomains_x);
    offsets[12] = offsetof(config_t, num_subdomains_y);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, mpi_config_t);
    MPI_Type_commit(mpi_config_t);
}

void mpi_shutdown(MPI_Datatype *mpi_config_t)
{
    MPI_Type_free(mpi_config_t);
    MPI_Finalize();
}

/* Sets the block type for the current processes
    ------------------------------
    |         |         |        |
    |    0    |    1    |    2   |
    |         |         |        |
    ------------------------------
    |         |         |        |
    |    3    |    4    |    5   |
    |         |         |        |
    ------------------------------
    |         |         |        |
    |    6    |    7    |    8   |
    |         |         |        |
    ------------------------------ */
int get_block_type(int rank, int num_subdomains_y, int num_subdomains_x)
{
    int x_loc, y_loc, block_type;

    x_loc = rank % num_subdomains_x;
    y_loc = (int)(rank / num_subdomains_y);
    block_type = 4;

    /* Checks if there is one subdomain */
    if ((num_subdomains_x == 1) && (num_subdomains_y == 1)) {
        block_type = 9;
        return block_type;
    }

    /* Corners */
    if ((x_loc == 0) && (y_loc == 0)) {
        block_type = 0;
    } else if ((x_loc == (num_subdomains_x - 1)) && (y_loc == 0)) {
        block_type = 2;
    } else if ((x_loc == 0) && (y_loc == (num_subdomains_y - 1))) {
        block_type = 6;
    } else if ((x_loc == (num_subdomains_x - 1)) && (y_loc == (num_subdomains_y - 1))) {
        block_type = 8;
    }

    /* Boundaries */
    else if (y_loc == 0) {
        block_type = 1;
    } else if (y_loc == (num_subdomains_y - 1)) {
        block_type = 7;
    } else if (x_loc == 0) {
        block_type = 3;
    } else if (x_loc == (num_subdomains_x - 1)) {
        block_type = 5;
    }

    return block_type;
}

/* Initializes the send the receive vectors */
void init_send_receive(send_vectors_t *send_vec, receive_vectors_t *rec_vec)
{
    /* Allocates memory for send and receive vectors */
    send_vec->send_vec_0 = malloc(dim.xdim * sizeof(double));
    send_vec->send_vec_1 = malloc(dim.ydim * sizeof(double));
    send_vec->send_vec_2 = malloc(dim.xdim * sizeof(double));
    send_vec->send_vec_3 = malloc(dim.ydim * sizeof(double));

    rec_vec->receive_vec_0 = malloc(dim.xdim * sizeof(double));
    rec_vec->receive_vec_1 = malloc(dim.ydim * sizeof(double));
    rec_vec->receive_vec_2 = malloc(dim.xdim * sizeof(double));
    rec_vec->receive_vec_3 = malloc(dim.ydim * sizeof(double));
}

/* Communication for type 0 block */
void comm_0(cell_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    printf("ranks: %d, %d, %d\n", rank, rank + 1, rank + dim.num_subdomains_x);

    /* Fills send vectors */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX((dim.ydim - 1), i)];
        send_vec->send_vec_2[i] = cur_cell->robin[2];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, (dim.xdim - 1))];
        send_vec->send_vec_1[i] = cur_cell->robin[1];
    }

    /* Sends send vectors */
    MPI_Send(send_vec->send_vec_1, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_2, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD);

    /* Receives rec vectors */
    MPI_Recv(rec_vec->receive_vec_3, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_0, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* Fills buffer region */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(dim.ydim, i)];
        cur_cell->robin[0] = rec_vec->receive_vec_0[i];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, dim.xdim)];
        cur_cell->robin[3] = rec_vec->receive_vec_3[i];
    }
}

/* Communication for type 1 block */
void comm_1(cell_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    /* Fills send vectors */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX((dim.ydim - 1), i)];
        send_vec->send_vec_2[i] = cur_cell->robin[2];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, (dim.xdim - 1))];
        send_vec->send_vec_1[i] = cur_cell->robin[1];

        cur_cell = &mesh[MESH_INDEX(i, 0)];
        send_vec->send_vec_3[i] = cur_cell->robin[3];
    }

    /* Sends send vectors */
    MPI_Send(send_vec->send_vec_1, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_2, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_3, dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);

    /* Receives rec vectors */
    MPI_Recv(rec_vec->receive_vec_3, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_0, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_1, dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* Fills buffer region */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(dim.ydim, i)];
        cur_cell->robin[0] = rec_vec->receive_vec_0[i];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, dim.xdim)];
        cur_cell->robin[3] = rec_vec->receive_vec_3[i];

        cur_cell = &mesh[MESH_INDEX(i, -1)];
        cur_cell->robin[1] = rec_vec->receive_vec_1[i];
    }
}

/* Communication for type 2 block */
void comm_2(cell_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    printf("ranks: %d, %d, %d\n", rank, rank - 1, rank + dim.num_subdomains_x);

    /* Fills send vectors */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX((dim.ydim - 1), i)];
        send_vec->send_vec_2[i] = cur_cell->robin[2];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, (dim.xdim - 1))];
        send_vec->send_vec_3[i] = cur_cell->robin[3];
    }

    /* Sends send vectors */
    MPI_Send(send_vec->send_vec_3, dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_2, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD);

    /* Receives rec vectors */
    MPI_Recv(rec_vec->receive_vec_1, dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_0, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* Fills buffer region */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(dim.ydim, i)];
        cur_cell->robin[0] = rec_vec->receive_vec_0[i];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, -1)];
        cur_cell->robin[1] = rec_vec->receive_vec_1[i];
    }
}

/* Communication for type 3 block */
void comm_3(cell_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    /* Fills send vectors */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX((dim.ydim - 1), i)];
        send_vec->send_vec_2[i] = cur_cell->robin[2];

        cur_cell = &mesh[MESH_INDEX(0, i)];
        send_vec->send_vec_0[i] = cur_cell->robin[0];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, (dim.xdim - 1))];
        send_vec->send_vec_1[i] = cur_cell->robin[1];
    }

    /* Sends send vectors */
    MPI_Send(send_vec->send_vec_1, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_2, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_0, dim.ydim, MPI_DOUBLE, rank - dim.num_subdomains_x, 0, MPI_COMM_WORLD);

    /* Receives rec vectors */
    MPI_Recv(rec_vec->receive_vec_3, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_0, dim.xdim, MPI_DOUBLE, rank + dim.num_subdomains_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_2, dim.ydim, MPI_DOUBLE, rank - dim.num_subdomains_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* Fills buffer region */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(dim.ydim, i)];
        cur_cell->robin[0] = rec_vec->receive_vec_0[i];

        cur_cell = &mesh[MESH_INDEX(-1, i)];
        cur_cell->robin[2] = rec_vec->receive_vec_2[i];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, dim.xdim)];
        cur_cell->robin[3] = rec_vec->receive_vec_3[i];
    }
}

/* Communication for type 6 block */
void comm_6(cell_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    printf("ranks: %d, %d, %d\n", rank, rank + 1, rank - dim.num_subdomains_x);

    /* Fills send vectors */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(0, i)];
        send_vec->send_vec_0[i] = cur_cell->robin[0];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, (dim.xdim - 1))];
        send_vec->send_vec_1[i] = cur_cell->robin[1];
    }

    /* Sends send vectors */
    MPI_Send(send_vec->send_vec_1, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_0, dim.ydim, MPI_DOUBLE, rank - dim.num_subdomains_x, 0, MPI_COMM_WORLD);

    /* Receives rec vectors */
    MPI_Recv(rec_vec->receive_vec_3, dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_2, dim.ydim, MPI_DOUBLE, rank - dim.num_subdomains_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* Fills buffer region */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(-1, i)];
        cur_cell->robin[2] = rec_vec->receive_vec_2[i];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, dim.xdim)];
        cur_cell->robin[3] = rec_vec->receive_vec_3[i];
    }
}

/* Communication for type 8 block */
void comm_8(cell_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    printf("ranks: %d, %d, %d\n", rank, rank - 1, rank - dim.num_subdomains_x);

    /* Fills send vectors */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(0, i)];
        send_vec->send_vec_0[i] = cur_cell->robin[0];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, 0)];
        send_vec->send_vec_3[i] = cur_cell->robin[3];
    }

    /* Sends send vectors */
    MPI_Send(send_vec->send_vec_3, dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
    MPI_Send(send_vec->send_vec_0, dim.xdim, MPI_DOUBLE, rank - dim.num_subdomains_x, 0, MPI_COMM_WORLD);

    /* Receives rec vectors */
    MPI_Recv(rec_vec->receive_vec_1, dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(rec_vec->receive_vec_2, dim.xdim, MPI_DOUBLE, rank - dim.num_subdomains_x, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* Fills buffer region */
    for (i = 0; i < dim.xdim; i++) {
        cur_cell = &mesh[MESH_INDEX(-1, i)];
        cur_cell->robin[2] = rec_vec->receive_vec_2[i];
    }

    for (i = 0; i < dim.ydim; i++) {
        cur_cell = &mesh[MESH_INDEX(i, -1)];
        cur_cell->robin[1] = rec_vec->receive_vec_1[i];
    }
}

/* Communication controller */
void comm(cell_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
            int block_type, int rank)
{
    switch (block_type) {
        case 0:
            comm_0(mesh, send_vec, rec_vec, rank);
            break;
        case 1:
            comm_1(mesh, send_vec, rec_vec, rank);
            break;
        case 2:
            comm_2(mesh, send_vec, rec_vec, rank);
            break;
        case 6:
            comm_6(mesh, send_vec, rec_vec, rank);
            break;
        case 8:
            comm_8(mesh, send_vec, rec_vec, rank);
            break;
    }
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi_util.h"

#define INDEX(y, x) (y * (mesh->dim.xdim + 2) + x)

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
int mpi_get_block_type(int rank, int num_subdomains_y, int num_subdomains_x)
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
void mpi_init_send_receive(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec)
{
    /* Allocates memory for send and receive vectors */
    send_vec->send_vec_0 = malloc(mesh->dim.xdim * sizeof(double));
    send_vec->send_vec_1 = malloc(mesh->dim.ydim * sizeof(double));
    send_vec->send_vec_2 = malloc(mesh->dim.xdim * sizeof(double));
    send_vec->send_vec_3 = malloc(mesh->dim.ydim * sizeof(double));

    rec_vec->receive_vec_0 = malloc(mesh->dim.xdim * sizeof(double));
    rec_vec->receive_vec_1 = malloc(mesh->dim.ydim * sizeof(double));
    rec_vec->receive_vec_2 = malloc(mesh->dim.xdim * sizeof(double));
    rec_vec->receive_vec_3 = malloc(mesh->dim.ydim * sizeof(double));
}

/* Sends Robin conditions to the right */
static void send_right(mesh_t *mesh, send_vectors_t *send_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), mesh->dim.xdim)];
        send_vec->send_vec_1[i] = cur_cell->robin[1];
    }

    MPI_Send(send_vec->send_vec_1, mesh->dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
}

/* Sends Robin conditions down */
static void send_down(mesh_t *mesh, send_vectors_t *send_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX(mesh->dim.ydim, (i + 1))];
        send_vec->send_vec_2[i] = cur_cell->robin[2];
    }

    MPI_Send(send_vec->send_vec_2, mesh->dim.xdim, MPI_DOUBLE, rank + mesh->dim.num_subdomains_x,
            0, MPI_COMM_WORLD);
}

/* Sends Robin conditions up */
static void send_up(mesh_t *mesh, send_vectors_t *send_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX(1, (i + 1))];
        send_vec->send_vec_0[i] = cur_cell->robin[0];
    }

    MPI_Send(send_vec->send_vec_0, mesh->dim.ydim, MPI_DOUBLE, rank - mesh->dim.num_subdomains_x,
        0, MPI_COMM_WORLD);
}

/* Sends Robin conditions left */
static void send_left(mesh_t *mesh, send_vectors_t *send_vec, int rank)
{
    int i;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), 1)];
        send_vec->send_vec_3[i] = cur_cell->robin[3];
    }

    MPI_Send(send_vec->send_vec_3, mesh->dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
}

/* Receive Robin conditions from right */
static void rec_right(mesh_t *mesh, receive_vectors_t *rec_vec, int rank)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_3, mesh->dim.ydim, MPI_DOUBLE, rank + 1, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), (mesh->dim.xdim + 1))];
        cur_cell->robin[3] = rec_vec->receive_vec_3[i];
    }
}

/* Receive Robin conditions from left */
static void rec_left(mesh_t *mesh, receive_vectors_t *rec_vec, int rank)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_1, mesh->dim.ydim, MPI_DOUBLE, rank - 1, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), 0)];
        cur_cell->robin[1] = rec_vec->receive_vec_1[i];
    }
}

/* Receive Robin conditions from down */
static void rec_down(mesh_t *mesh, receive_vectors_t *rec_vec, int rank)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_0, mesh->dim.xdim, MPI_DOUBLE, rank + mesh->dim.num_subdomains_x,
        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX((mesh->dim.ydim + 1), (i + 1))];
        cur_cell->robin[0] = rec_vec->receive_vec_0[i];
    }
}

/* Receive Robin conditions from down */
static void rec_up(mesh_t *mesh, receive_vectors_t *rec_vec, int rank)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_2, mesh->dim.xdim, MPI_DOUBLE, rank - mesh->dim.num_subdomains_x,
        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX(0, (i + 1))];
        cur_cell->robin[2] = rec_vec->receive_vec_2[i];
    }
}

/* Communication for type 0 block */
static void comm_0(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_right(mesh, send_vec, rank);
    send_down(mesh, send_vec, rank);

    rec_right(mesh, rec_vec, rank);
    rec_down(mesh, rec_vec, rank);
}

/* Communication for type 1 block */
static void comm_1(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_left(mesh, send_vec, rank);
    send_right(mesh, send_vec, rank);
    send_down(mesh, send_vec, rank);

    rec_left(mesh, rec_vec, rank);
    rec_right(mesh, rec_vec, rank);
    rec_down(mesh, rec_vec, rank);
}

/* Communication for type 2 block */
static void comm_2(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_left(mesh, send_vec, rank);
    send_down(mesh, send_vec, rank);

    rec_left(mesh, rec_vec, rank);
    rec_down(mesh, rec_vec, rank);
}

/* Communication for type 3 block */
static void comm_3(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_right(mesh, send_vec, rank);
    send_up(mesh, send_vec, rank);
    send_down(mesh, send_vec, rank);

    rec_right(mesh, rec_vec, rank);
    rec_up(mesh, rec_vec, rank);
    rec_down(mesh, rec_vec, rank);
}

/* Communication for type 4 block */
static void comm_4(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_right(mesh, send_vec, rank);
    send_up(mesh, send_vec, rank);
    send_down(mesh, send_vec, rank);
    send_left(mesh, send_vec, rank);

    rec_right(mesh, rec_vec, rank);
    rec_up(mesh, rec_vec, rank);
    rec_down(mesh, rec_vec, rank);
    rec_left(mesh, rec_vec, rank);
}

/* Communication for type 5 block */
static void comm_5(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_up(mesh, send_vec, rank);
    send_down(mesh, send_vec, rank);
    send_left(mesh, send_vec, rank);

    rec_up(mesh, rec_vec, rank);
    rec_down(mesh, rec_vec, rank);
    rec_left(mesh, rec_vec, rank);
}

/* Communication for type 6 block */
static void comm_6(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_right(mesh, send_vec, rank);
    send_up(mesh, send_vec, rank);

    rec_right(mesh, rec_vec, rank);
    rec_up(mesh, rec_vec, rank);
}

/* Communication for type 7 block */
static void comm_7(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_right(mesh, send_vec, rank);
    send_up(mesh, send_vec, rank);
    send_left(mesh, send_vec, rank);

    rec_right(mesh, rec_vec, rank);
    rec_up(mesh, rec_vec, rank);
    rec_left(mesh, rec_vec, rank);
}

/* Communication for type 8 block */
static void comm_8(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec, int rank)
{
    send_left(mesh, send_vec, rank);
    send_up(mesh, send_vec, rank);

    rec_left(mesh, rec_vec, rank);
    rec_up(mesh, rec_vec, rank);
}

/* Communication controller */
void mpi_comm(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
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
        case 3:
            comm_3(mesh, send_vec, rec_vec, rank);
            break;
        case 4:
            comm_4(mesh, send_vec, rec_vec, rank);
            break;
        case 5:
            comm_5(mesh, send_vec, rec_vec, rank);
            break;
        case 6:
            comm_6(mesh, send_vec, rec_vec, rank);
            break;
        case 7:
            comm_7(mesh, send_vec, rec_vec, rank);
            break;
        case 8:
            comm_8(mesh, send_vec, rec_vec, rank);
            break;
    }
}

/* Recombines Pressure Output */
void mpi_recombine_pressure(mesh_t *mesh, int rank)
{
    int i, j;
    cell_t *cur_cell;
}

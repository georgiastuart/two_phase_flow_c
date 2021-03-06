#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include "mpi_util.h"
#include "mesh.h"

#define INDEX(y, x) (y * (mesh->dim.xdim + 2) + x)
#define INDEX_NO_MESH(y, x, xdim) (y * (xdim) + x)
#define NDIMS 2
#define CONFIG_LEN 28
#define STR_LEN 100
#define ARRAY_SPLIT 9999;

void mpi_setup(int *argc, char ***argv, int *rank, int *size, MPI_Datatype *mpi_config_t)
{
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, size);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);

    /* Create type for config struct */
    const int nitems = CONFIG_LEN;
    int blocklengths[CONFIG_LEN] = {1, 1, 1, 1, STR_LEN, STR_LEN, 1, 1, 1, 1, 1, 1, 1,
                            STR_LEN, STR_LEN, STR_LEN, 1, 1, 1, 1, 1, 1, STR_LEN,
                            1, 1, 1, STR_LEN, STR_LEN};
    MPI_Datatype types[CONFIG_LEN] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE,
                                MPI_CHAR, MPI_CHAR,
                                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                MPI_INT, MPI_INT, MPI_INT, MPI_CHAR, MPI_CHAR, MPI_CHAR,
                                MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                                MPI_DOUBLE, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_CHAR,
                                MPI_CHAR};
    MPI_Aint offsets[CONFIG_LEN];

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
    offsets[13] = offsetof(config_t, pressure_out);
    offsets[14] = offsetof(config_t, velocity_y_out);
    offsets[15] = offsetof(config_t, velocity_x_out);
    offsets[16] = offsetof(config_t, porosity);
    offsets[17] = offsetof(config_t, sat_rel_o);
    offsets[18] = offsetof(config_t, sat_rel_w);
    offsets[19] = offsetof(config_t, visc_o);
    offsets[20] = offsetof(config_t, visc_w);
    offsets[21] = offsetof(config_t, eta);
    offsets[22] = offsetof(config_t, saturation_out);
    offsets[23] = offsetof(config_t, time_steps);
    offsets[24] = offsetof(config_t, dt);
    offsets[25] = offsetof(config_t, linearity);
    offsets[26] = offsetof(config_t, sat_file);
    offsets[27] = offsetof(config_t, prod_well_out);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, mpi_config_t);
    MPI_Type_commit(mpi_config_t);
}

/* For transmitting permeability and source information to the different processes */
/* Mode is 0 for permeability, 1 for source, 2 for saturation */
void mpi_setup_parameters(config_t *config, int mode, int size, int is_master, double **out_param)
{
    int i, j, k, xdim_per_block, ydim_per_block, x_block_loc, y_block_loc;
    int transmit_len;
    double *transmit, *full_param;
    // MPI_Request request = MPI_REQUEST_NULL;

    xdim_per_block = config->xdim / config->num_subdomains_x;
    ydim_per_block = config->ydim / config->num_subdomains_y;

    transmit_len = (xdim_per_block + 2) * (ydim_per_block + 2);

    /* Allocate memory for transmit vector and parameter vector */
    transmit = malloc(transmit_len * sizeof(double));
    *out_param = malloc(transmit_len * sizeof(double));
    full_param = malloc((config->xdim + 2) * (config->ydim + 2) * sizeof(double));

    if (is_master) {
        if (mode == 0) {
            full_param = read_file_pad(config->perm_file, config->ydim, config->xdim);
        } else if (mode == 1) {
            full_param = read_file_pad(config->src_file, config->ydim, config->xdim);
        } else if (mode == 2) {
            full_param = read_file_pad(config->sat_file, config->ydim, config->xdim);
        } else {
            printf("Invalid Mode\n");
            exit(1);
        }

        for (k = 1; k < size; k++) {
            x_block_loc = k % config->num_subdomains_x;
            y_block_loc = (int)(k / config->num_subdomains_y);

            for (i = 0; i < (ydim_per_block + 2); i++) {
                for (j = 0; j < (xdim_per_block + 2); j++) {
                    transmit[INDEX_NO_MESH(i, j, xdim_per_block + 2)] =
                                full_param[INDEX_NO_MESH((i + y_block_loc * ydim_per_block),
                                (j + x_block_loc * xdim_per_block), (config->xdim + 2))];
                }
            }
            MPI_Ssend(transmit, transmit_len, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(*out_param, transmit_len, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    if (is_master) {
        k = 0;
        x_block_loc = k % config->num_subdomains_x;
        y_block_loc = (int)(k / config->num_subdomains_y);

        for (i = 0; i < (ydim_per_block + 2); i++) {
            for (j = 0; j < (xdim_per_block + 2); j++) {
                (*out_param)[INDEX_NO_MESH(i, j, xdim_per_block + 2)] =
                            full_param[INDEX_NO_MESH((i + y_block_loc * ydim_per_block),
                            (j + x_block_loc * xdim_per_block), (config->xdim + 2))];
            }
        }
    }
    free(transmit);
    free(full_param);
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

/* Sends values to the right */
/* Mode is 0 for Robin, 1 for Saturation */
static void send_right(mesh_t *mesh, send_vectors_t *send_vec, int rank, int mode)
{
    int i;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), mesh->dim.xdim)];
        if (mode) {
            send_vec->send_vec_1[i] = cur_cell->saturation;
        } else {
            send_vec->send_vec_1[i] = cur_cell->robin[1];
        }
    }

    MPI_Send(send_vec->send_vec_1, mesh->dim.ydim, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
}

/* Sends values down */
/* Mode is 0 for Robin, 1 for Saturation */
static void send_down(mesh_t *mesh, send_vectors_t *send_vec, int rank, int mode)
{
    cell_t *cur_cell;

    for (int i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX(mesh->dim.ydim, (i + 1))];
        if (mode) {
            send_vec->send_vec_2[i] = cur_cell->saturation;
        } else {
            send_vec->send_vec_2[i] = cur_cell->robin[2];
        }
    }

    MPI_Send(send_vec->send_vec_2, mesh->dim.xdim, MPI_DOUBLE, rank + mesh->dim.num_subdomains_x,
            0, MPI_COMM_WORLD);
}

/* Sends Robin conditions up */
static void send_up(mesh_t *mesh, send_vectors_t *send_vec, int rank, int mode)
{
    int i;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX(1, (i + 1))];
        if (mode) {
            send_vec->send_vec_0[i] = cur_cell->saturation;
        } else {
            send_vec->send_vec_0[i] = cur_cell->robin[0];
        }
    }

    MPI_Send(send_vec->send_vec_0, mesh->dim.ydim, MPI_DOUBLE, rank - mesh->dim.num_subdomains_x,
        0, MPI_COMM_WORLD);
}

/* Sends Robin conditions left */
static void send_left(mesh_t *mesh, send_vectors_t *send_vec, int rank, int mode)
{
    int i;
    cell_t *cur_cell;

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), 1)];
        if (mode) {
            send_vec->send_vec_3[i] = cur_cell->saturation;
        } else {
            send_vec->send_vec_3[i] = cur_cell->robin[3];
        }
    }

    MPI_Send(send_vec->send_vec_3, mesh->dim.ydim, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
}

/* Receive Robin conditions from right */
static void rec_right(mesh_t *mesh, receive_vectors_t *rec_vec, int rank, int mode)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_3, mesh->dim.ydim, MPI_DOUBLE, rank + 1, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), (mesh->dim.xdim + 1))];
        if (mode) {
            cur_cell->saturation = rec_vec->receive_vec_3[i];
        } else {
            cur_cell->robin[3] = rec_vec->receive_vec_3[i];
        }
    }
}

/* Receive Robin conditions from left */
static void rec_left(mesh_t *mesh, receive_vectors_t *rec_vec, int rank, int mode)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_1, mesh->dim.ydim, MPI_DOUBLE, rank - 1, 0,
        MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.ydim; i++) {
        cur_cell = &mesh->cell[INDEX((i + 1), 0)];
        if (mode) {
            cur_cell->saturation = rec_vec->receive_vec_1[i];
        } else {
            cur_cell->robin[1] = rec_vec->receive_vec_1[i];
        }
    }
}

/* Receive Robin conditions from down */
static void rec_down(mesh_t *mesh, receive_vectors_t *rec_vec, int rank, int mode)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_0, mesh->dim.xdim, MPI_DOUBLE, rank + mesh->dim.num_subdomains_x,
        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX((mesh->dim.ydim + 1), (i + 1))];
        if (mode) {
            cur_cell->saturation = rec_vec->receive_vec_0[i];
        } else {
            cur_cell->robin[0] = rec_vec->receive_vec_0[i];
        }
    }
}

/* Receive Robin conditions from down */
static void rec_up(mesh_t *mesh, receive_vectors_t *rec_vec, int rank, int mode)
{
    cell_t *cur_cell;
    int i;

    MPI_Recv(rec_vec->receive_vec_2, mesh->dim.xdim, MPI_DOUBLE, rank - mesh->dim.num_subdomains_x,
        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (i = 0; i < mesh->dim.xdim; i++) {
        cur_cell = &mesh->cell[INDEX(0, (i + 1))];
        if (mode) {
            cur_cell->saturation = rec_vec->receive_vec_2[i];
        } else {
            cur_cell->robin[2] = rec_vec->receive_vec_2[i];
        }
    }
}

/* Communication for type 0 block */
static void comm_0(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_right(mesh, send_vec, rank, mode);
    send_down(mesh, send_vec, rank, mode);

    rec_right(mesh, rec_vec, rank, mode);
    rec_down(mesh, rec_vec, rank, mode);
}

/* Communication for type 1 block */
static void comm_1(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_left(mesh, send_vec, rank, mode);
    send_right(mesh, send_vec, rank, mode);
    send_down(mesh, send_vec, rank, mode);

    rec_left(mesh, rec_vec, rank, mode);
    rec_right(mesh, rec_vec, rank, mode);
    rec_down(mesh, rec_vec, rank, mode);
}

/* Communication for type 2 block */
static void comm_2(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_left(mesh, send_vec, rank, mode);
    send_down(mesh, send_vec, rank, mode);

    rec_left(mesh, rec_vec, rank, mode);
    rec_down(mesh, rec_vec, rank, mode);
}

/* Communication for type 3 block */
static void comm_3(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_right(mesh, send_vec, rank, mode);
    send_up(mesh, send_vec, rank, mode);
    send_down(mesh, send_vec, rank, mode);

    rec_right(mesh, rec_vec, rank, mode);
    rec_up(mesh, rec_vec, rank, mode);
    rec_down(mesh, rec_vec, rank, mode);
}

/* Communication for type 4 block */
static void comm_4(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_right(mesh, send_vec, rank, mode);
    send_up(mesh, send_vec, rank, mode);
    send_down(mesh, send_vec, rank, mode);
    send_left(mesh, send_vec, rank, mode);

    rec_right(mesh, rec_vec, rank, mode);
    rec_up(mesh, rec_vec, rank, mode);
    rec_down(mesh, rec_vec, rank, mode);
    rec_left(mesh, rec_vec, rank, mode);
}

/* Communication for type 5 block */
static void comm_5(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_up(mesh, send_vec, rank, mode);
    send_down(mesh, send_vec, rank, mode);
    send_left(mesh, send_vec, rank, mode);

    rec_up(mesh, rec_vec, rank, mode);
    rec_down(mesh, rec_vec, rank, mode);
    rec_left(mesh, rec_vec, rank, mode);
}

/* Communication for type 6 block */
static void comm_6(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_right(mesh, send_vec, rank, mode);
    send_up(mesh, send_vec, rank, mode);

    rec_right(mesh, rec_vec, rank, mode);
    rec_up(mesh, rec_vec, rank, mode);
}

/* Communication for type 7 block */
static void comm_7(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_right(mesh, send_vec, rank, mode);
    send_up(mesh, send_vec, rank, mode);
    send_left(mesh, send_vec, rank, mode);

    rec_right(mesh, rec_vec, rank, mode);
    rec_up(mesh, rec_vec, rank, mode);
    rec_left(mesh, rec_vec, rank, mode);
}

/* Communication for type 8 block */
static void comm_8(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
                    int rank, int mode)
{
    send_left(mesh, send_vec, rank, mode);
    send_up(mesh, send_vec, rank, mode);

    rec_left(mesh, rec_vec, rank, mode);
    rec_up(mesh, rec_vec, rank, mode);
}

/* Communication controller */
/* Mode 0 for Robin, 1 for Saturation */
void mpi_comm(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
            int block_type, int rank, int mode)
{
    switch (block_type) {
        case 0:
            comm_0(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 1:
            comm_1(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 2:
            comm_2(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 3:
            comm_3(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 4:
            comm_4(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 5:
            comm_5(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 6:
            comm_6(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 7:
            comm_7(mesh, send_vec, rec_vec, rank, mode);
            break;
        case 8:
            comm_8(mesh, send_vec, rec_vec, rank, mode);
            break;
    }
}

void write_data(mesh_t *mesh, config_t *config, int size, int rank, const char *mode)
{
    int gsizes[NDIMS], distribs[NDIMS], dargs[NDIMS], psizes[NDIMS], file_type_size;
    MPI_Datatype file_type;
    MPI_Aint file_type_extent;
    MPI_File fh;
    int write_buffer_size;
    double *write_buffer;
    int i, j, file_open_error, file_write_error;
    MPI_Status status;
    char name[100];

    /* Allocates memory for write buffer */
    write_buffer  = malloc(mesh->dim.xdim * mesh->dim.ydim * sizeof(double));

    /* Transfers data from mesh to write buffer */
    if (!strcmp(mode, "pressure")) {
        strcpy(name, config->pressure_out);
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                write_buffer[MESH_INDEX_NO_PAD(i, j)] = mesh->cell[MESH_INDEX(i, j)].pressure;
            }
        }
    } else if (!strcmp(mode, "velocity_y")) {
        strcpy(name, config->velocity_y_out);
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                write_buffer[MESH_INDEX_NO_PAD(i, j)] = mesh->cell[MESH_INDEX(i, j)].velocity_y;
            }
        }
    } else if (!strcmp(mode, "velocity_x")) {
        strcpy(name, config->velocity_x_out);
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                write_buffer[MESH_INDEX_NO_PAD(i, j)] = mesh->cell[MESH_INDEX(i, j)].velocity_x;
            }
        }
    } else if (!strcmp(mode, "saturation")) {
        strcpy(name, config->saturation_out);
        for (i = 0; i < mesh->dim.ydim; i++) {
            for (j = 0; j < mesh->dim.xdim; j++) {
                write_buffer[MESH_INDEX_NO_PAD(i, j)] = mesh->cell[MESH_INDEX(i, j)].saturation;
            }
        }
    } else {
        printf("Invalid print mode\n");
        exit(1);
    }

    /* Sets up MPI darray */
    gsizes[0] = config->ydim;
    gsizes[1] = config->xdim;
    distribs[0] = MPI_DISTRIBUTE_BLOCK;
    distribs[1] = MPI_DISTRIBUTE_BLOCK;
    dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
    psizes[0] = config->num_subdomains_y;
    psizes[1] = config->num_subdomains_x;

    MPI_Type_create_darray(size, rank, NDIMS, gsizes, distribs, dargs, psizes,
        MPI_ORDER_C, MPI_DOUBLE, &file_type);
    MPI_Type_commit(&file_type);

    MPI_Type_extent(file_type, &file_type_extent);
    MPI_Type_size(file_type, &file_type_size);

    write_buffer_size = mesh->dim.xdim * mesh->dim.ydim;

    MPI_Barrier(MPI_COMM_WORLD);

    /* Deletes old file */
    MPI_File_delete(name, MPI_INFO_NULL);

    /* Opens write out file */
    file_open_error = MPI_File_open(MPI_COMM_WORLD, name,
				    MPI_MODE_CREATE | MPI_MODE_WRONLY,
				    MPI_INFO_NULL, &fh);

    if (file_open_error != MPI_SUCCESS) {
        printf("Rank %d failed to open file", rank);
        MPI_Abort(MPI_COMM_WORLD, file_open_error);
        exit(1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Sets file view */
    MPI_File_set_view(fh, 0, MPI_DOUBLE, file_type, "native", MPI_INFO_NULL);

    /* Writes file out */
    file_write_error = MPI_File_write_all(fh, write_buffer, write_buffer_size,
                        MPI_DOUBLE, &status);

    if (file_write_error != MPI_SUCCESS) {
        printf("Rank %d file write error\n", rank);

        MPI_File_close(&fh);
        free(write_buffer);
        if (rank == 0)
            MPI_File_delete(name, MPI_INFO_NULL);
        exit(1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Finished writing %s file\n", mode);
    }

    MPI_File_close(&fh);
    free(write_buffer);
}

/* Writes production wells */
void write_prod_well_data(production_wells_t *wells, config_t *config, int size, int rank)
{
    MPI_File fh;
    double *buffer;
    int buffer_size, num_doubles;
    int offset, send_offset, index_base;

    /* Allocates memory for buffer */
    num_doubles = config->time_steps + 3;
    buffer_size = wells->num_wells * (num_doubles) * sizeof(double);
    buffer = malloc(buffer_size);

    /* Deletes old file */
    MPI_File_delete(config->prod_well_out, MPI_INFO_NULL);

    /* Opens file */
    int file_open_error = MPI_File_open(MPI_COMM_WORLD, config->prod_well_out,
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    if (file_open_error != MPI_SUCCESS) {
        printf("Rank %d file read error\n", rank);
        exit(1);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Sets up offsets */
    offset = 0;

    if (rank != 0) {
        MPI_Recv(&offset, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    send_offset = offset + buffer_size;

    if (rank != (size - 1)) {
        MPI_Send(&send_offset, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }

    /* Reads data into buffer */
    for (int i = 0; i < wells->num_wells; i++) {
        index_base = i * num_doubles;
        buffer[0 + index_base] = (double) wells->wells[i].y_pos;
        buffer[1 + index_base] = (double) wells->wells[i].x_pos;

        for (int j = 0; j < config->time_steps; j++) {
            buffer[2 + j + index_base] = wells->wells[i].oil_sat_recording[j];
            // printf("Buffer: %e\n", buffer[2 + j + index_base]);
            // printf("Oil sat rec: %e\n", wells->wells[i].oil_sat_recording[j]);
        }

        buffer[2 + config->time_steps + index_base] = ARRAY_SPLIT;
    }

    // printf("num doubles %d\n", num_doubles);
    /* Writes data to file at specified offset */
    if (wells->num_wells > 0) {
        int file_write_error = MPI_File_write_at(fh, offset, buffer, num_doubles,
                                MPI_DOUBLE, MPI_STATUS_IGNORE);

        if (file_write_error != MPI_SUCCESS) {
            printf("Rank %d failed to write file\n", rank);
            exit(1);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Finished writing production well file\n");
    }

    MPI_File_close(&fh);

    free(buffer);
}

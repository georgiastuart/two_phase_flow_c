#ifndef H_MPI_UTIL
#define H_MPI_UTIL

#include "util.h"
#include "mesh.h"

typedef struct send_vectors
{
    double *send_vec_0, *send_vec_1, *send_vec_2, *send_vec_3;
} send_vectors_t;

typedef struct receive_vectors
{
    double *receive_vec_0, *receive_vec_1, *receive_vec_2, *receive_vec_3;
} receive_vectors_t;

void mpi_setup(int *argc, char ***argv, int *rank, int *size, MPI_Datatype *mpi_config_t);
void mpi_shutdown(MPI_Datatype *mpi_config_t);
void mpi_init_send_receive(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec);
void mpi_comm(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
            int block_type, int rank);
int mpi_get_block_type(int rank, int num_subdomains_y, int num_subdomains_x);


#endif  /* H_MPI_UTIL */

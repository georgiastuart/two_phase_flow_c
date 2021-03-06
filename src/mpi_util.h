#ifndef H_MPI_UTIL
#define H_MPI_UTIL

#include "util.h"

struct mesh;
typedef struct mesh mesh_t;

struct production_wells;
typedef struct production_wells production_wells_t;

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
void mpi_setup_parameters(config_t *config, int mode, int size, int is_master, double **param);
void mpi_init_send_receive(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec);
void mpi_comm(mesh_t *mesh, send_vectors_t *send_vec, receive_vectors_t *rec_vec,
            int block_type, int rank, int mode);
int mpi_get_block_type(int rank, int num_subdomains_y, int num_subdomains_x);
void write_data(mesh_t *mesh, config_t *config, int size, int rank, const char *mode);
void write_prod_well_data(production_wells_t *wells, config_t *config, int size, int rank);


#endif  /* H_MPI_UTIL */

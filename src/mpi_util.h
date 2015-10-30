#ifndef H_MPI_UTIL
#define H_MPI_UTIL

#include "mpi.h"
#include "util.h"

void mpi_setup(int *argc, char ***argv, int *rank, int *size, MPI_Datatype *mpi_config_t);
void mpi_shutdown(MPI_Datatype *mpi_config_t);


#endif  /* H_MPI_UTIL */

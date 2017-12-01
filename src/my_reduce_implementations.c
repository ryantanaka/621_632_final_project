#include <mpi.h>
#include <math.h>

// Num_ranks must be power of 2!
int Reduce_binomial(void *sendbuf_notype, void* recvbuf_notype, int m,
MPI_Datatype mpi_datatype, MPI_Op mpi_op, int root, MPI_Comm mpi_comm) {
  int my_rank, num_ranks, num_phases, i, k;
  MPI_Status status;

  int *sendbuf, *recvbuf, *tempbuf;

  sendbuf = (int *) sendbuf_notype;
  recvbuf = (int *) recvbuf_notype;

  MPI_Comm_rank(mpi_comm, &my_rank);
  MPI_Comm_size(mpi_comm, &num_ranks);

  tempbuf = (int*)malloc(m * sizeof(int));

  if (my_rank == 0) {
    for (k = 0; k < m; k++) {
      recvbuf[k] = sendbuf[k];
    }
  }
  else {
    for (k = 0; k < m; k++) {
      tempbuf[k] = sendbuf[k];
    }
  }

  num_phases = (int)log2(num_ranks);

  for (i = 0; i < num_phases; i++) {
    if (my_rank == 0) {
      MPI_Recv(tempbuf, m, mpi_datatype, my_rank + pow(2, i), 0, mpi_comm, &status);

      for (k = 0; k < m; k++) {
        recvbuf[k] += tempbuf[k];
      }
    }
    else if (my_rank % (int)(pow(2, i+1)) == 0) {
      MPI_Recv(tempbuf, m, mpi_datatype, my_rank + pow(2, i), 0, mpi_comm, &status);

      for (k = 0; k < m; k++) {
        sendbuf[k] += tempbuf[k];
      }
    }
    else if (my_rank % (int)pow(2, i) == 0) {
      MPI_Send(sendbuf, m, mpi_datatype, my_rank - pow(2, i), 0, mpi_comm);
    }
  }

  free(tempbuf);

  MPI_Barrier(mpi_comm);
}

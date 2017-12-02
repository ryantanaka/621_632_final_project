#include <mpi.h>
#include <math.h>
#include <stdio.h>

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
}

// segment size s must be (1 <= s <= m)
int Reduce_pipeline(void *sendbuf_notype, void* recvbuf_notype, int s, int m,
MPI_Datatype mpi_datatype, MPI_Op mpi_op, int root, MPI_Comm mpi_comm) {
  int my_rank, num_ranks;
  int i, j;
  int s2, s3, snext, q;
  int send_to, receive_from, start_rank;
  MPI_Status status;

  int *sendbuf, *recvbuf, *tempbuf;

  sendbuf = (int *) sendbuf_notype;
  recvbuf = (int *) recvbuf_notype;

  MPI_Comm_rank(mpi_comm, &my_rank);
  MPI_Comm_size(mpi_comm, &num_ranks);

  send_to = (my_rank + num_ranks - 1) % num_ranks;
  receive_from = (my_rank + num_ranks + 1) % num_ranks;
  start_rank = (root + num_ranks - 1) % num_ranks;

  tempbuf = (int*)malloc(s * sizeof(int));

  if (my_rank == root) {
    for (i = 0; i < m; i++) {
      recvbuf[i] = sendbuf[i];
    }
  }
  else {
    for(i = 0; i < s; i++) {
      tempbuf[i] = sendbuf[i];
    }
  }

  q = m/s; // q is the number of chunks
  if ((m % s) != 0) {
    s2 = m % s;
    q++;
  }
  else {
    s2 = s;
  }

  // DEBUGGING
  if (my_rank == root) {
    printf("segment size s = %d, segment size 2 s2 = %d, num chunks q = %d\n", s, s2, q);
  }

  // im only sending chunks
  if (my_rank == start_rank) {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Send(tempbuf, s3, mpi_datatype, send_to, 0, mpi_comm);
      if (i < q - 1) { //get the next chunk ready to send if there are chunks left to send
        snext = (i == q - 2) ? s2 : s;
        for(j = 0; j < snext; j++) {
          tempbuf[j] = sendbuf[(i + 1) * s + j];
        }
      }
    }
  }
  // im a rank in between start_rank and root
  // receiving, reducing, then sending
  else if (my_rank != root) {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Recv(tempbuf, s3, mpi_datatype, receive_from, 0, mpi_comm, &status);
      for (j = 0; j < s3; j++) {
        tempbuf[j] += sendbuf[i * s + j];
      }
      MPI_Send(tempbuf, s3, mpi_datatype, send_to, 0, mpi_comm);
    }
  }
  // im the root
  // im only receiving then reducing chunks
  else {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Recv(tempbuf, s3, mpi_datatype, receive_from, 0 , mpi_comm, &status);
      for(j = 0; j < s3; j++) {
        recvbuf[i*s + j] += tempbuf[j];
      }
    }
  }

  free(tempbuf);
}

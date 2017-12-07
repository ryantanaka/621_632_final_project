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

  tempbuf = (int*)malloc(m*sizeof(int));

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
int Reduce_pipeline(void *sendbuf_notype, void* recvbuf_notype, int m,
MPI_Datatype mpi_datatype, MPI_Op mpi_op, int root, MPI_Comm mpi_comm) {
  int my_rank, num_ranks;
  int i, j;
  int s2, s3, snext, q;
  int s = m/2;
  int send_to, receive_from, start_rank;
  MPI_Status status;

  int *sendbuf, *recvbuf, *tempbuf, *tempbuf2;

  sendbuf = (int *) sendbuf_notype;
  recvbuf = (int *) recvbuf_notype;

  MPI_Comm_rank(mpi_comm, &my_rank);
  MPI_Comm_size(mpi_comm, &num_ranks);

  send_to = (my_rank + num_ranks - 1) % num_ranks;
  receive_from = (my_rank + num_ranks + 1) % num_ranks;
  start_rank = (root + num_ranks - 1) % num_ranks;

  if(my_rank != root) {
    tempbuf = (int*)malloc(s*sizeof(int));
  }
  tempbuf2 = (int*)malloc(s*sizeof(int));

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
  /*
  if (my_rank == root) {
    printf("segment size s = %d, segment size 2 s2 = %d, num chunks q = %d\n", s, s2, q);
  }
  */

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
  else if (my_rank != start_rank && my_rank != root) {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Recv(tempbuf2, s3, mpi_datatype, receive_from, 0, mpi_comm, &status);
      for (j = 0; j < s3; j++) {
        tempbuf[j] += tempbuf2[j];
      }
      /* FOR DEBUGGING
      if(my_rank == 1) {
        for(j=0; j<m; j++) {
          printf(" %d", tempbuf[j]);
        }
      }
      */
      MPI_Send(tempbuf, s3, mpi_datatype, send_to, 0, mpi_comm);
      if (i < q - 1) { //get the next chunk ready to send if there are chunks left to send
        snext = (i == q - 2) ? s2 : s;
        for(j = 0; j < snext; j++) {
          tempbuf[j] = sendbuf[(i + 1) * s + j];
        }
      }
    }
  }
  else {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Recv(tempbuf2, s3, mpi_datatype, receive_from, 0 , mpi_comm, &status);
      for(j = 0; j < s3; j++) {
        recvbuf[i*s + j] += tempbuf2[j];
      }
    }
  }

  if(my_rank != root) {
    free(tempbuf);
  }
  free(tempbuf2);
}

int Reduce_binary(void *sendbuf_notype, void* recvbuf_notype, int m,
MPI_Datatype mpi_datatype, MPI_Op mpi_op, int root, MPI_Comm mpi_comm) {
  int my_rank, num_ranks;
  int i,j;
  int s2, s3, snext, q;
  int s = m/2;
  int send_to, receive_from;
  int left_child, right_child, parent, is_leaf;
  int *tree_layout;
  int my_tree_index;
  int *sendbuf, *recvbuf, *tempbuf, *tempbuf2;
  MPI_Status status;
  MPI_Request request1 = MPI_REQUEST_NULL;

  sendbuf = (int *) sendbuf_notype;
  recvbuf = (int *) recvbuf_notype;

  MPI_Comm_rank(mpi_comm, &my_rank);
  MPI_Comm_size(mpi_comm, &num_ranks);

  if(my_rank != root) {
    tempbuf = (int*)malloc(s*sizeof(int));
  }
  tempbuf2 = (int*)malloc(s*sizeof(int));

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

  // START computing node info

  // build level order tree based on whatever the root is (root will be at index 0)
  tree_layout = (int *)malloc(num_ranks * sizeof(int));
  for (i = 0; i < num_ranks; i++) {
    tree_layout[i] = (i + root) % num_ranks;
    if (tree_layout[i] == my_rank) {
      my_tree_index = i;
    }
  }

  // figure out parent and children
  if(my_rank != root) {
    parent = tree_layout[(int)ceil((double)my_tree_index / (double)2) - 1];
  }

  left_child = 2*my_tree_index + 1;
  right_child = 2*my_tree_index + 2;

  if (left_child >= num_ranks) {
    left_child = -1;
  }
  else {
    left_child = tree_layout[left_child];
  }

  if (right_child >= num_ranks) {
    right_child = -1;
  }
  else {
    right_child = tree_layout[right_child];
  }

  if (left_child == -1 && right_child == -1) {
    is_leaf = 1;
  }
  else {
    is_leaf = 0;
  }

  free(tree_layout);
  // END computing node info

  if (my_rank == root) {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Recv(tempbuf2, s3, mpi_datatype, left_child, 0 , mpi_comm, &status);
      for(j = 0; j < s3; j++) {
        recvbuf[i * s + j] += tempbuf2[j];
      }

      if(right_child != -1) {
        MPI_Recv(tempbuf2, s3, mpi_datatype, right_child, 0 , mpi_comm, &status);
        for(j = 0; j < s3; j++) {
          recvbuf[i * s + j] += tempbuf2[j];
        }
      }
    }
  }
  else if(is_leaf) {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Send(tempbuf, s3, mpi_datatype, parent, 0, mpi_comm);
      if (i < q - 1) { //get the next chunk ready to send if there are chunks left to send
        snext = (i == q - 2) ? s2 : s;
        for(j = 0; j < snext; j++) {
          tempbuf[j] = sendbuf[(i + 1) * s + j];
        }
      }
    }
  }
  else {
    for(i = 0; i < q; i++) {
      s3 = (i == q - 1) ? s2 : s;

      MPI_Recv(tempbuf2, s3, mpi_datatype, left_child, 0 , mpi_comm, &status);
      for(j = 0; j < s3; j++) {
        tempbuf[j] += tempbuf2[j];
      }

      if(right_child != -1) {
        MPI_Recv(tempbuf2, s3, mpi_datatype, right_child, 0 , mpi_comm, &status);
        for(j = 0; j < s3; j++) {
          tempbuf[j] += tempbuf2[j];
        }
      }

      MPI_Send(tempbuf, s3, mpi_datatype, parent, 0, mpi_comm);
      if (i < q - 1) { //get the next chunk ready to send if there are chunks left to send
        snext = (i == q - 2) ? s2 : s;
        for(j = 0; j < snext; j++) {
          tempbuf[j] = sendbuf[(i + 1) * s + j];
        }
      }
    }
  }

  if(my_rank != root) {
    free(tempbuf);
  }
  free(tempbuf2);
}

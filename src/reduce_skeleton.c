#include <unistd.h>
#include <stdio.h>
#include <mpi.h>

#define NUM_INTS 10000000

int Reduce_greedy( void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int Reduce_binomial( void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);

int main(int argc, char **argv) {
  // IMPORTANT: each process has its own i, my_rank, nprocs, x
  int i, my_rank, num_procs;
  //int local_sum[5] = {1,2,3,4,5};
  int *local_sum, *global_sum;

  local_sum = malloc(sizeof(int) * NUM_INTS);
  global_sum = malloc(sizeof(int) * NUM_INTS);
  for (int k = 0; k < NUM_INTS; k++) {
    local_sum[k] = 1;
  }

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  // Start the timer
  double start_time;
  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) {
    start_time = MPI_Wtime();
  }

  MPI_Reduce(local_sum, global_sum, NUM_INTS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  //Reduce_binomial(local_sum, global_sum, NUM_INTS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  //Reduce_greedy(local_sum, global_sum, NUM_INTS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  // End timer and print result
  MPI_Barrier (MPI_COMM_WORLD);
  if (my_rank == 0) {
    printf("REDUCE TIME: %.3lf seconds\n", MPI_Wtime() - start_time);

    /*
    for (int j = 0; j < NUM_INTS; j++) {
      printf("%d ", global_sum[j]);
    }
    printf("\n");
    */
  }

  MPI_Finalize();
}

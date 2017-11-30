#include <unistd.h>
#include <stdio.h>
#include <mpi.h>

int Reduce_greedy( void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);
int Reduce_binomial( void *, void *, int, MPI_Datatype, MPI_Op, int, MPI_Comm);

int main(int argc, char **argv) {
  int i, my_rank, num_procs;
  int *local_sum, *global_sum;
  char *reduce_implementation_name;

  local_sum = malloc(sizeof(int) * NUM_INTS);
  global_sum = malloc(sizeof(int) * NUM_INTS);
  for (int k = 0; k < NUM_INTS; k++) {
    local_sum[k] = 1;
  }

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  reduce_implementation_name = argv[1];

  // Start the timer
  double start_time;
  MPI_Barrier(MPI_COMM_WORLD);
  if (my_rank == 0) {
    start_time = MPI_Wtime();
  }

  if (strcmp(reduce_implementation_name, "smpi_binomial_reduce") == 0) {
    MPI_Reduce(local_sum, global_sum, NUM_INTS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  else if (strcmp(reduce_implementation_name, "binomial_reduce") == 0) {
    Reduce_binomial(local_sum, global_sum, NUM_INTS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  else if (strcmp(reduce_implementation_name, "greedy_reduce") == 0){
    Reduce_greedy(local_sum, global_sum, NUM_INTS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  // End timer and print result
  MPI_Barrier (MPI_COMM_WORLD);
  if (my_rank == 0) {
    printf("%s | NUM_INTS: %d | %.3lf seconds\n", reduce_implementation_name, NUM_INTS, MPI_Wtime() - start_time);

    /*
    for (int j = 0; j < NUM_INTS; j++) {
      printf("%d ", global_sum[j]);
    }
    printf("\n");
    */
  }

  MPI_Finalize();
}

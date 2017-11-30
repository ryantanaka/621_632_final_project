#include <mpi.h>

// segment size
int global_s = 50000;

int Reduce_greedy( void *sendbuf_notype, void *recvbuf_notype, int m,
MPI_Datatype mpi_datatype, MPI_Op mpi_op, int root, MPI_Comm mpi_comm) {

  int pool_size; // num_ranks
  int my_rank;
  int i; // hcurrent segment index
  int j; // buffer index
  int q; // num segments
  int z; // offset to determine where to send / receive (ie. send to (my_rank - z))
  int s2, s3; // hold segment sizes
  int qstart, qend; // segment we start working from, segment we end on
  int snext; // size of next segment
  int s=global_s; // segment size
  MPI_Status status;
  int *tempbuf, *tempbuf2; // hold send / recieve values before sending or receiving
  int *sendbuf, *recvbuf; // pointers to local sum and global sum
  int *hist; // array that says where we have sent a specific segment

  sendbuf = (int *) sendbuf_notype;
  recvbuf = (int *) recvbuf_notype;

  if (root != 0 ) {
    MPI_Abort( mpi_comm, 512);
  }

  MPI_Comm_size(mpi_comm, &pool_size);
  MPI_Comm_rank(mpi_comm, &my_rank);

  if(my_rank != 0 ) {
    tempbuf = (int*)malloc(s*sizeof(int));
  }
  tempbuf2 = (int*)malloc(s*sizeof(int));

  if ( my_rank == 0) {
    for(i=0;i<m;i++) {
      recvbuf[i] = sendbuf[i];
    }
  }
  // for everyone else put sendbuf into tempbuf
  else {
    for(i=0;i<s;i++) {
      tempbuf[i] = sendbuf[i];
    }
  }

  q = m/s; // we have q segements
  if( (m % s) != 0 ) { // if total message size isnt nicely divdied by segment size
    s2 = m % s; // we have a segment2 = the remainder
    q++; // because of the remainder we have 1 more full segment
  }
  else {
    s2 = s;  // last segment is same as earlier segments (they are all equal size)
  }

  // table to keep track of each segment
  hist = (int *)malloc(q*sizeof(int));
  for(i = 0; i < q; i++) {
    hist[i] = pool_size; // set each segment in the table
  }                      // to the number of procs involved

  qstart = 0; // start of segment(s) to work on (this is index of the segment in segments, not specific ints)
  qend = 0; // end of segment(s) to work on

  // while the last segment hasn't been touched by p-1 procs??
  while( hist[q-1] > 1 ) {
    // we can start a new segment if the current one has been touched by two or more procs compared to the next??
    if (qstart != q-1 && hist[qstart+1]-hist[qstart] > 1) { // HERE: do oldest first!
      qstart++;
    }

    for(i = qstart; i >= qend; i--) { // for the segments we are currently working on
      z = (i == 0) ? hist[0]/2 : (hist[i]-hist[i-1])/2; // z will determine where we send / receive message
      s3 = ( i == q-1 ) ? s2 : s; // s3 will be the segment size you send (s for regular or s2 for the remainder segment)

      if( my_rank < hist[i] && my_rank >= hist[i] - z ) {
        MPI_Send(tempbuf, s3, mpi_datatype, my_rank-z, 512, mpi_comm );
        if( i < q-1 ) {
          snext = (i == q-2) ? s2 : s; // switch to handle when the last segment is of different size
          for(j=0;j<snext;j++) {
            tempbuf[j] = sendbuf[(i+1)*s + j]; // put the next segment into tempbuf so it can be sent
          }                                    // if we are working on any segment before the last
        }
      }

      if( my_rank < hist[i] - z && my_rank >= hist[i] - 2*z ) { // you are receiving if you are on the lower half of procs
        MPI_Recv( tempbuf2, s3, mpi_datatype, my_rank+z, 512, mpi_comm, &status );
        if ( my_rank == 0 ) {
          for( j = 0;j < s3; j++ ) {
            recvbuf[i*s + j] += tempbuf2[j]; // reduce
          }
        }
        else {
          for( j = 0; j < s3; j++) {
            tempbuf[j] += tempbuf2[j]; // reduce
          }
        }
      }
      hist[i] -= z; // decrement hist by z
      if (hist[i] == 1) { // when we are done with a segment 
        hist[i] = 0;
        qend++; // move qend index up one because segment is done
      }
    }
  }

  free( hist );

  if(my_rank != 0) {
    free( tempbuf );
  }

  free( tempbuf2 );
  return 0;
}


#include <mpi.h>
#include <iostream>
using namespace std;
#define  ARRAYSIZE	200000000

double  data[ARRAYSIZE];

int main (int argc, char *argv[])
{
  int numtasks,
      rank,
      rc,
      offset,
      tag1 = 2,
      tag2 = 1,
      source;

double rank_sum, sum;
double update(int myoffset, int chunk, int myid);

MPI_Status status;

double start_time, end_time;
start_time = MPI_Wtime();

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);

cout << "MPI task - " << rank << " of " << numtasks << " tasks" << endl;

int chunksize = (ARRAYSIZE / numtasks);
int leftover = (ARRAYSIZE % numtasks);

if (rank == 0){
  sum = 0;
  offset = chunksize + leftover;

  for(int i=0; i<ARRAYSIZE; i++) 
    {
    data[i] =  (i + 1.0);
    }

  for (int dest = 1; dest < numtasks; dest++) 
  {
    MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
    MPI_Send(&data[offset], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
    offset = offset + chunksize;
  }

  offset = 0;
  rank_sum = update(offset, chunksize+leftover, rank);

  for (int i=1; i< numtasks; i++) {
    source = i;
    MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
    MPI_Recv(&data[offset], chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);
    }

    MPI_Reduce(&rank_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    cout << "Final result = " << sum << endl;
    end_time   = MPI_Wtime();
    cout << "Time =" << end_time - start_time << endl;
  } 

if (rank > 0) {

  int source = 0;
  MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
  MPI_Recv(&data[offset], chunksize, MPI_DOUBLE, source, tag2, MPI_COMM_WORLD, &status);

  double start_rank_time = MPI_Wtime();
  rank_sum = update(offset, chunksize, rank);
  double end_rank_time = MPI_Wtime();
  cout << "rank " << rank <<" time = " << end_rank_time - start_rank_time << endl << endl;

  int dest = 0;
  MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
  MPI_Send(&data[offset], chunksize, MPI_DOUBLE, 0, tag2, MPI_COMM_WORLD);

  // MPI_REDUCE(sendbuf, recvbuf, count, datatype, op, root, comm)
  MPI_Reduce(&rank_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  } 

  MPI_Finalize();

}  

double update(int offset, int chunk, int rank) {
  int i; 
  double sum =  0.0;
  for(int i =offset;i< offset+chunk;i++)
  {
      sum += 1/data[i];
  }
  cout << "rank = " << rank << " sum = " << sum <<endl;
  return(sum);
  }

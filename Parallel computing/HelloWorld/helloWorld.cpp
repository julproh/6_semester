#include <mpi.h>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
  
    // Get the number of processes ssociated with the communicator
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the calling process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    cout << "Hello world from process - "<< processor_name <<  " with rank: " << world_rank << " out of " << world_size << " processors\n"<< endl;

    // Finalize: Any resources allocated for MPI can be freed
    MPI_Finalize();
    return 0;
}
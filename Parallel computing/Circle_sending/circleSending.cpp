#include <iostream>
#include <mpi.h>

using namespace std;

int main(int argc, char* argv[]) {

	int process;
	int rank;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	// Функция определения номера процесса MPI_Comm_rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Функция определения числа процессов в области связи MPI_Comm_size
	MPI_Comm_size(MPI_COMM_WORLD, &process);

	if (rank == 0) {

		int data = 0;

		// int MPI_Send(void* message, int count,MPI_Datatype datatype, int dest, int tag,MPI_Comm comm)
		MPI_Send(&data, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);

		cout << "From " << rank << ": send " << data << endl;

		MPI_Recv(&data, 1, MPI_INT, process - 1, 1, MPI_COMM_WORLD, &status);

		cout << "From  " << rank << ": receive " << data << endl;
	}
	else {

		int data;

		MPI_Recv(&data, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, &status);

		data +=rank;

		if (rank == process - 1) { MPI_Send(&data, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); }
		else { MPI_Send(&data, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD); }

		cout << "From " << rank << ": receive " << data - rank << ", send " << data << endl;
	}

	MPI_Finalize();

	return 0;
}
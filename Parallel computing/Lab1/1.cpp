#include<stdio.h>
#include<iostream>
#include<mpi.h>
#include <iomanip>
#include<vector>

using namespace std;

#define SIZE_X 10
#define SIZE_T 10


int M = SIZE_X; // по х
int K = SIZE_T; // по t

//будем рассматривать линейный случай
double tau = 0.01;
double h = 0.01;

double phi(double x) {
	return x;
}

double ksi(double t) {
	return t;
}

double f(double t, double x) {
	return 0;
}

//Решение крестом
double solve(double left, double right, double bottom, double f) {
	return bottom + 2 * tau * (f + (left - right) / (2 * h));
}

//Решение явной центральной трехточечной схемой на нижней границе
double solve_bottom(double left, double right, double f) {
	return 0.5 * (left + right) + tau * (f + (left - right) / (2 * h));
}

//Решение уголком на правой границе
double solve_right(double left, double central, double f) {
	return central + tau * (f + (left - central) / h);
}

//отправкка данных
void data_send(double** data, MPI_Status& status, int ProcNum, int ProcRank, int N_PerProcess, int k);

//прием данных
pair<double, double> data_recive(double** data, MPI_Status& status, int ProcNum, int ProcRank, int N_PerProcess, int k);

//выделение памяти под матрицу
template <typename T>
T** new_matrix(int height, int width) {
	T** matrix = new T * [height];
	for (int i = 0; i < height; ++i) {
		matrix[i] = new T[width];
	}
	return matrix;
}

//удаление матрицы
template <typename T>
void delete_matrix(T** matrix, int height) {
	for (int i = 0; i < height; ++i) {
		delete[] matrix[i];
	}
	delete[] matrix;
}

int main(int argc, char* argv[]) {

	int ProcNum;
	int ProcRank;

	MPI_Status status;

	double time_d = 0;
	double time_end = 0;
	
	int N_PerProcess, N_Add;
	int FirstX;

	double time_start = MPI_Wtime();
	// Инициализируем среду MPI
	MPI_Init(&argc, &argv); 
	// Опереляем ранг текущего процесса
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	// Определяем общее количество процессов
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);

	if (ProcRank == 0) {
		cout << "From process " << ProcRank << ": calculations begin" << endl;
	}

	//определяем количество точек для каждого процесса
	N_PerProcess = M / ProcNum;
	N_Add = M % ProcNum;

	if (ProcRank < N_Add) {
		N_PerProcess++;
		//определяем номер элеманта, с которого начнутся вычисления этого процесса
		FirstX = N_PerProcess * ProcRank;
	}
	else {
		FirstX = N_Add + N_PerProcess * ProcRank;
	}

	double** data = new_matrix<double>(K, N_PerProcess);

	// Процесс вычисляет свою часть значений в 1ой строке
	for (int n = 0; n < N_PerProcess; ++n) data[0][n] = phi((n + FirstX) * h);

	// Обмен данными для полноценного заполнения первой строки
	data_send(data, status, ProcNum, ProcRank, N_PerProcess, 0);

	// Получение  значения 1ой строки
	pair<double, double> data_recv = data_recive(data, status, ProcNum, ProcRank, N_PerProcess, 0);
	
	// Заполнение середины 2ой строки (также по частям)
	for (int n = 1; n < N_PerProcess - 1; ++n) data[1][n] = solve_bottom(data[0][n - 1], data[0][n + 1], f(tau, (FirstX + n) * h));

	// Определение первого значения (своей группы) в 2ой строке
	if (ProcRank == 0)
		data[1][0] = ksi(tau);
	else				
		data[1][0] = solve_bottom(data_recv.first, data[0][1], f(tau, FirstX * h));

	// Определение последнего значения (своей группы) 2ой строки
	if (ProcRank == ProcNum - 1)	
		data[1][N_PerProcess - 1] = solve_right(data[0][N_PerProcess - 2], data[0][N_PerProcess - 1], f(tau, h * (FirstX + N_PerProcess - 1)));
	else
		data[1][N_PerProcess - 1] = solve_bottom(data[0][N_PerProcess - 2], data_recv.second, f(tau, h * (FirstX + N_PerProcess - 1)));

	// Отступили от края - можно считать для всех остальных
	for (int k = 2; k < K; ++k) {
		
		// Пересылка данных
		data_send(data, status, ProcNum, ProcRank, N_PerProcess, k - 1);

		// Значения в середине
		for (int n = 1; n < N_PerProcess - 1; ++n) data[k][n] = solve(data[k - 1][n - 1], data[k - 2][n], data[k - 1][n + 1], f(k * tau, (FirstX + n) * h));

		// Получение данных
		pair<double, double> data_recv = data_recive(data, status, ProcNum, ProcRank, N_PerProcess, k - 1);

		// Граничное условие слева
		if (ProcRank == 0) data[k][0] = ksi(k * tau);
		else data[k][0] = solve(data_recv.first, data[k - 2][0], data[k - 1][1], f(k * tau, FirstX));

		// Граничное условие справа 
		if (ProcRank == ProcNum - 1) data[k][N_PerProcess - 1] = solve_right(data[k - 1][N_PerProcess - 2], data[k - 1][N_PerProcess - 1], f(k * tau, (FirstX + N_PerProcess - 1) * h));
		else data[k][N_PerProcess - 1] = solve(data[k - 1][N_PerProcess - 2], data[k - 2][N_PerProcess - 1], data_recv.second, f(k * tau, (FirstX + N_PerProcess - 1) * h));
	}

	//Сбор посчитанных данных в процесс 0
	if (ProcRank != 0) {
		double* reviving_data = new double[N_PerProcess * K];

		// Подготовка массива для данных
		for (int y = 0; y < K; ++y)
			for (int x = 0; x < N_PerProcess; ++x) reviving_data[y * N_PerProcess + x] = data[y][x];

		// Отправление количество элеменнтов, обрабатываемых процессом
		MPI_Send(&N_PerProcess, sizeof(N_PerProcess), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
		// Отправление номера элемента, с которого начинал
		MPI_Send(&FirstX, sizeof(FirstX), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
		// Отправление данных
		MPI_Send(reviving_data, N_PerProcess * K * sizeof(double), MPI_BYTE, 0, 1, MPI_COMM_WORLD);

		cout << "From process " << ProcRank << ": data sended to process 0 succesfully " << endl;

		delete[] reviving_data;
	}
	else {
		double** finale_data = new_matrix<double>(K, M);
		for (int y = 0; y < K; ++y)
			for (int x = 0; x < N_PerProcess; ++x)
				finale_data[y][x] = data[y][x];

		// Принимаем данных
		for (int i = 1; i < ProcNum; ++i) {

			// Колличество обрабатываемых столбцов у процесса, который пересылает данные
			int rec_N_PerProcess;
			MPI_Recv(&rec_N_PerProcess, sizeof(rec_N_PerProcess), MPI_BYTE, i, 1, MPI_COMM_WORLD, &status);

			// Номер элемента, с которого начанался рассчет
			int proc_FirstX;
			MPI_Recv(&proc_FirstX, sizeof(proc_FirstX), MPI_BYTE, i, 1, MPI_COMM_WORLD, &status);

			// Значения
			double* reciving_data = new double[rec_N_PerProcess * K];
			MPI_Recv(reciving_data, rec_N_PerProcess * K * sizeof(double), MPI_BYTE, i, 1, MPI_COMM_WORLD, &status);

			cout << "From process " << ProcRank << ": data recived from process  " << i << " succesfully" << endl;

			for (int iy = 0; iy < K; ++iy)
				for (int ix = 0; ix < rec_N_PerProcess; ++ix)
					finale_data[iy][ix + proc_FirstX] = reciving_data[iy * rec_N_PerProcess + ix];
		}

		// Затраченное время
		time_end = MPI_Wtime();
		time_d = time_end - time_start;

		cout << "From process " << ProcRank << ": calculations completed. Time " << time_d << " seconds" << endl;

		cout << endl <<  "Final data: "  << endl << endl;

		cout << fixed;
		cout.precision(2);

		// Вывод рассчитанных данных
		for (int iy = 0; iy < K; ++iy) {
			for (int ix = 0; ix < M; ++ix) {
				cout << setw(10) << finale_data[iy][ix] << " ";
			}
			cout << endl;
		}

		delete_matrix<double>(finale_data, K);
	}

	delete_matrix<double>(data, K);

	MPI_Finalize();

	return 0;
}

//отправкка данных
void data_send(double** data, MPI_Status& status, int ProcNum, int ProcRank, int N_PerProcess, int k) {
	
	double send_right, send_left;

	if (ProcRank == 0) {
		send_right = data[k][N_PerProcess - 1];
		MPI_Send(&send_right, sizeof(send_right), MPI_BYTE, 1, 1, MPI_COMM_WORLD);
	}
	else if (ProcRank == ProcNum - 1) {
		send_left = data[k][0];
		MPI_Send(&send_left, sizeof(send_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD);
	}
	else {
		send_right = data[k][N_PerProcess - 1];
		MPI_Send(&send_right, sizeof(send_right), MPI_BYTE, ProcRank + 1, 1, MPI_COMM_WORLD);
		send_left = data[k][0];
		MPI_Send(&send_left, sizeof(send_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD);
	}
}

//прием данных
pair<double, double> data_recive(double** data, MPI_Status& status, int ProcNum, int ProcRank, int N_PerProcess, int k) {
	double recv_right, recv_left;

	if (ProcRank == 0) {
		MPI_Recv(&recv_right, sizeof(recv_right), MPI_BYTE, 1, 1, MPI_COMM_WORLD, &status);
		recv_left = 0;
	}
	else if (ProcRank == ProcNum - 1) {
		MPI_Recv(&recv_left, sizeof(recv_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD, &status);
		recv_right = 0;
	}
	else {
		MPI_Recv(&recv_right, sizeof(recv_right), MPI_BYTE, ProcRank + 1, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&recv_left, sizeof(recv_left), MPI_BYTE, ProcRank - 1, 1, MPI_COMM_WORLD, &status);
	}

	return pair<double, double>(recv_left, recv_right);
}
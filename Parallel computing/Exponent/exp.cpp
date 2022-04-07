#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <mpi.h>
using namespace std;

// Вычисление факториала
long long fact (int n) {
    long long res = 1;
    for (int i = 1; i <= n; i++) {
        res *= i;
    }
    return res;
}

int main(int argc, char* argv[]) {
    
    int N = atoi(argv[1]); // Из строки argv[1] находим число N
    MPI_Init(&argc, &argv);
    
    double start = MPI_Wtime();

    int proc_ttl, proc_num;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_ttl);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_num);
    double total_sum = 0, part_sum = 0;
    int k = N / proc_ttl;
    int left, right;
    double proc_start = MPI_Wtime();

    // Определение границ 
    if (proc_num == 0) {
        left = 1;
        right = N - (proc_ttl - 1) * k;
    }
    else {
        left = N - proc_num * k + 1;
        right = N - (proc_num - 1) * k;
    }

    long long curr = fact(left);
    for (int i = left; i < right; i++) {
        // Вычисление знаменателя i-го члена ряда
        if (i != left) {
            curr *= i;
        }
        // Суммирование с оставшимся рядом
        part_sum += ((double) 1) / curr;
    }

    // Сбор подсумм в proc_num = 0
    MPI_Reduce(&part_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    double proc_finish = MPI_Wtime();
  
    cout << "Time for rank " << proc_num << " is " << proc_finish - proc_start << endl;
    double finish = MPI_Wtime();
   
    if (proc_num == 0) {
        // добавление 0-го члена (+ 1)
        cout << "exp = " << setprecision(12) << 1 + total_sum << endl;
        cout << "Total time: " << finish - start << endl;
    }
    MPI_Finalize();
    return 0;
}
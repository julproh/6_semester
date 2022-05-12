/* Условие
С помощью pthread необходимо сделать программу, которая параллельно вычислит определенный интеграл sin(1/x)
 в пределе от некого положительного действительного числа до некоторого большего положительного действительного числа. 
1) Необходимо сбалансировать решение - время выполнения задачи на всех нитях должно быть одинаковым.
2) Шаги интегрирования должны быть динамическими.
*/

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <pthread.h>
#include <ctime>
using namespace std;

//Определяем количество потоков
#define NUM_THREADS 4

//Точное значение интеграла
#define _VALUE 0.504066 

//Инициализация мьютекса
pthread_mutex_t mutex_ = PTHREAD_MUTEX_INITIALIZER;

double VALUE = 0;

typedef struct thread_data {
	double start = 0.0;
	double end = 0.0;
	double eps = 0.0;
	double step = 0.0;
	int id = 0;
} thread_data;

double func(double x) {
	return sin(1.0/x);
}

//Реализация формулы трапеции 
double trapezoidal(double n, double a, double b);
double count_integral(double start, double end, double eps, double step);

//Расспределение нагрузки
double* balanced(double a, double b, int thread_num);
double* balanced_inside(double a, double b, double eps, double* step, int thread_num);

//Функция для каждого потока
void* count_thread(void* args);

int main(int argc, char** argv) {
	
	double eps = 0.02;
	double a = 0.001;
	double b = 1.0;

	//Рабиваем отрезок на промежутки
	double* points = balanced(a, b, NUM_THREADS);

	for (int i = 0; i < NUM_THREADS; i++)
	{
		double step = 0;
		//Внетренние промежутки делим на части
		double* points_inside = balanced_inside(points[i], points[i + 1], eps, &step,  NUM_THREADS);
	
		pthread_t threads[NUM_THREADS];
		thread_data td[NUM_THREADS];

		for(int j = 0; j < NUM_THREADS; j++) {
                    
			td[j].start = points_inside[j];
			td[j].end = points_inside[j + 1];
			td[j].eps = eps;
			td[j].step = step;
			td[j].id = j;
			
		//Cоздание потока и передача ему массива
		if((pthread_create(&threads[j], NULL, count_thread, (void*)&td[j]))!=0)
		{
			cout << "error: pthread_create" << endl;
			return EXIT_FAILURE;
		}
		}

		//Ожидание завершения работы потоков
		for(int j = 0; j < NUM_THREADS; j++) 
			pthread_join(threads[j],NULL);
		free(points_inside);
	}

	free(points);

	cout << "Значение интергара = " << VALUE << endl;
	cout << "Ошибка от точного значения = " << fabs(_VALUE - VALUE) << endl;
	pthread_exit(NULL);
	return 0;
}

double trapezoidal(double n, double a, double b)
{
    double h = (b-a)/n;
    double s = func(a)+func(b);

    for (int i = 1; i < n; i++)
        s += 2*func(a+i*h);

    return (h/2)*s;
}

double count_integral(double start, double end, double eps, double step) {
	int n = (end - start)/step;
	return trapezoidal(2*n, start, end);
}

double* balanced(double a, double b, int thread_num) {

	double a_inv = 1 / a;
	double b_inv = 1 / b;
	double linv = a_inv - b_inv;
	double part_linv = linv / thread_num;

	double* points = (double*) calloc (thread_num+1, sizeof(double));

	cout << "Опорные точки: " ;

	for(int i = 0; i <= thread_num; i++)
	{
		points[i] = 1 / (a_inv - part_linv * i);
		cout << points[i] << " ";
	}
	cout << endl;

	return points;
}

double* balanced_inside(double a, double b, double eps, double* step, int thread_num) {

	(*step) = (b-a)/10;
	while(fabs(func(a)-func(a + (*step))) >= eps)
	{
		(*step) /= 10;
	}
	double* points = (double*) calloc (thread_num+1, sizeof(double));

	cout << "Точки: ";
	for(int i = 0; i <= thread_num; i++)
	{
		points[i] = a + i * (b - a)/thread_num;
		cout << points[i] << " ";
	}
	cout << endl;

	return points;
}

void* count_thread(void* args) {

	time_t time = clock();
	thread_data* data = (thread_data*)args;

	double integral = count_integral(data -> start, data -> end, data -> eps / NUM_THREADS / NUM_THREADS , data -> step);
	VALUE += integral;
	
    //Временно оганичиваем доступ для других потоков
    pthread_mutex_lock(&mutex_);
	cout << "Затраченное " <<  pthread_self() << " время : "  << (clock() - time) << " ms " 
			 << " на [ " <<(data -> start) <<  " ; " << data -> end << "] рассчитанный интеграл: " << integral << endl; 
			 // cout  << "Количетсво точек: " << (data -> end-data -> start)/(data -> step) ;  
	pthread_mutex_unlock(&mutex_);
}
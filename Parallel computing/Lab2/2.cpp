/* Условие
С помощью pthread необходимо сделать программу, которая параллельно вычислит определенный интеграл sin(1/x)
 в пределе от некого положительного действительного числа до некоторого большего положительного действительного числа. 
1) Необходимо сбалансировать решение - время выполнения задачи на всех нитях должно быть одинаковым.
2) Шаги интегрирования должны быть динамическими.
*/

#include <pthread.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#define NUM_THREADS 8

typedef struct _thread_data_t {
	int tid;
	int stuff;
	int div_num;
} thread_data_t;

double shared_summ;
pthread_mutex_t lock_x;


//globals
double initial_limit_a = 0.001;
double initial_limit_b = 1.0;
double accuracy = 0.000000001;



void *thr_func(void *arg);
double function(double x);


int main(int argc, char **argv)
{
	double error = (1.0 * 1e+12 + 2.0 * 1e+9) / 24 * 319;
	int divide_num = (int)((1 - 1.0 /M_PI) / pow(accuracy / error, 1.0/3)) + 1;
	//std::cout <<error<< " "<< (int)((1 - 1.0 /M_PI) / pow(accuracy / error, 1.0/3)) << " "<<pow((accuracy / error), (1.0/3))<<  std::endl;
	pthread_t thr[NUM_THREADS];
	int i, rc;
	thread_data_t thr_data[NUM_THREADS];
	shared_summ = 0;
	pthread_mutex_init(&lock_x, NULL);

	for(i = 0; i < NUM_THREADS; i++)
	{
		thr_data[i].tid = i;
		thr_data[i].stuff = i;
		thr_data[i].div_num = divide_num;
		if((rc = pthread_create(&thr[i], NULL, thr_func, &thr_data[i])))
		{
			printf("error: pthread_create, rc: %d\n", rc);
			return EXIT_FAILURE;
		}
 	}

	for(i = 0; i< NUM_THREADS; ++i)
	{
		pthread_join(thr[i], NULL);
	}

	std::cout << "Integral of sin(1/x) from 0.001 to 1 with accuracy " << accuracy << " is equal ";
	std::cout.precision(10);
	std::cout << shared_summ << std::endl;

	return EXIT_SUCCESS;
}

void *thr_func(void *arg)
{
	thread_data_t *data = (thread_data_t *)arg;
	double counter = data -> stuff + 1;
	double limit_a = 1.0 / M_PI / counter;
	double limit_b;
	int divide_num = data -> div_num;
	double summ = 0;

	while(limit_a >= 1.0 / M_PI / 319)
	{
		if (limit_a < initial_limit_a)
				limit_a = initial_limit_a;
		if (counter == 1)
			limit_b = initial_limit_b;
		else
			limit_b =  1.0 / M_PI / (counter - 1);

		//if(data->tid == 6)
		//	std::cout <<"["<< limit_a << ','<< limit_b << "] " << divide_num<< std::endl;

		for(int i = 0; i < divide_num; i++)
			summ += (limit_b - limit_a) / divide_num * function((limit_b - limit_a) / divide_num * i + limit_a);

		counter += NUM_THREADS;
		limit_a = 1.0 / M_PI / counter;
		limit_b =  1.0 / M_PI / (counter - 1);
	}

	pthread_mutex_lock(&lock_x);
	shared_summ += summ;
	pthread_mutex_unlock(&lock_x);
	pthread_exit(NULL);
}


double function(double x)
{
	return sin(1.0 / x);
}

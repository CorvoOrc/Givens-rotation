/*Steshenko Alexander*/

/*
This is my implementation Givens rotation method
used parallel technology MPI
http://en.wikipedia.org/wiki/Givens_rotation
*/

// input: dimension
// output: correctness of decision and work time

#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
//#include <iostream>
//#include <conio.h>
#include <time.h>
#include <mpi.h>

const double EPS = 0.00001;

double a[100000000];
double ans[10000];

int get_deep(int n) {
	int i = 0, j = 0;
	int deep = 0;

	int k = 2;
	for (j = k; k <= n * 2; j += k, ++deep) {
		if (j > n) {
			k *= 2, j = 0, --deep;
			continue;
		}

		if (i >= n - 1) {
			k *= 2, j = 0;
		}
	}

	return deep;
}

void init(int &n, int &rank) {
	//printf("rank=%d n=%d", rank, n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n + 1; ++j) {
			if (i == j) a[i*n + j] = 2.1 * (n - 1);
			else a[i*n + j] = 1.0;
		}

		ans[i] = a[i*n + n + 1] / a[i*n + i];
	}
}

bool check(int &n) {
	bool verdict = true;

	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum = 0;
		for (int j = 0; j < n; ++j) sum += a[i*n + j] * ans[j];
		//std::cout << sum << " " << b[i] << std::endl;
		if (abs(sum - a[i*n + n + 1]) >= EPS) {
			verdict = false;
			//printf("i=%d sum=%f b[i]=%f abs=%f\n", i, sum, b[i], abs(sum - b[i]));
			//std::cout << i << ". " << sum << " " << b[i] << " " << abs(sum - b[i]) << std::endl;
		}
	}

	return verdict;
}

int main(int argc, char* argv[])
{
	int n = atoi(argv[1]),
		my_rank = 9,
		n_threads = 9,
		i, j,
		first,
		last,
		comm,
		k = 0;
	double sum,
		start_time,
		finish_time;

	int *a_count = (int *)malloc(sizeof(int) * n),
		*a_number = (int *)malloc(sizeof(int) * n),
		*a_count2 = (int *)malloc(sizeof(int) * n),
		*a_number2 = (int *)malloc(sizeof(int) * n);

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &n_threads);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	comm = MPI_COMM_WORLD;

	if (my_rank == 0) {
		start_time = MPI_Wtime();
		init(n, my_rank);
	}

	first = my_rank*n / n_threads;
	last = (my_rank + 1)*n / n_threads - 1;

	for (i = 0; i < n_threads; ++i) {
		a_count[i] = (i + 1) * n / n_threads - i*n / n_threads;
		a_number[i] = i*n / n_threads;
	}

	for (i = 0; i < n_threads; ++i) {
		a_count2[i] = a_count[my_rank];
		a_number2[i] = a_number[my_rank];
	}

	MPI_Bcast(a, n*(n+1), MPI_DOUBLE_PRECISION, 0, comm);
	MPI_Bcast(ans, n, MPI_DOUBLE_PRECISION, 0, comm);

	// main loop
	do {
		int p = 2;
		i = k;
		int d = get_deep(n - k - 1); //(n - k - 1) % 2==0 ? n - k - 1 : n - k;
		int d_i = 0;

		for (i = first; n - k >= p / 2 && i <= last; i += p, ++d_i) {
			if (i >= n - 1) {
				p *= 2, i = k, --d_i;
				continue;
			}

			double dist = pow(a[k*n + k] * a[k*n + k] + a[(i - k / 2 - 1)*n + k] * a[(i - k / 2 - 1)*n + k], 0.5);
			double c = a[k*n + k] / dist,
				   s = a[(i - k / 2 - 1)*n + k] / dist;
			for (j = k; j <= n + 1; ++j) {
				double  old_value_a_ij = a[(i - k / 2 - 1)*n + j],
					    old_value_a_kj = a[k*n + j];

				a[k*n + j] = c * old_value_a_kj + s * old_value_a_ij;
				a[(i - k / 2 - 1)*n + j] = c * old_value_a_ij - s * old_value_a_kj;
			}

			if (i >= n - 2)
				p *= 2, i = k;
		}
		MPI_Barrier(comm);

		for (i = first; i <= last; ++i)
			MPI_Alltoallv(&a[i], a_count2, a_number2, MPI_DOUBLE_PRECISION, &a[i], a_count, a_number, MPI_DOUBLE_PRECISION, comm);

		++k;
	} while (k < n);

	MPI_Barrier(comm);

	if (my_rank == 0) {
		for (i = n; i > 0; --i) {
			sum = 0;
			for (j = n; j > i; --j)
				sum += ans[j] * a[i*n + j];

			ans[i] = (a[i*n + n + 1] - sum) / a[i*n + i]; 
		}

		finish_time = MPI_Wtime();
		printf("Ax=b? ");
		if (check(n) == 1)
			printf("Yes\n");
		else
			printf("No\n");
		printf("Time = %f\n", finish_time - start_time);
	}

	MPI_Finalize();

	//getch();
	return 0;
}

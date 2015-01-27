/* Steshenko Alexander */

/*
This is my implementation Givens rotation method
used parallel technology OpenMP
http://en.wikipedia.org/wiki/Givens_rotation
*/

// input: dimension, number threads
// output: correctness of decision and work time

#include <stdlib.h>
#include <iostream>
//#include <conio.h>
#include <math.h>
#include <omp.h>

using namespace std;

const double EPS = 0.00001;

void print_result(double *(&out), int n) {
	for (int i = 1; i < n + 1; ++i)
		std::cout << "x[" << i << "] = " << out[i] << std::endl;
}

int get_deep(int n) {
	int i = 0, j = 0;
	int deep = 0;

	int k = 2;
	for (j = k; k <= n*2; j += k, ++deep) {
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

void init(double *a, double *x, int &n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n+1; ++j) {
			if (i == j) a[i*n + j] = rand() / 100;
			else a[i*n + j] = rand() / 100;
		}

		x[i] = a[i*n + n + 1] / a[i*n + i];
		//x_new[i] = 0.0;
	}
}

bool check(double *a, double *x, int &n) {
	bool verdict = true;

	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum = 0;
		for (int j = 0; j < n; ++j) sum += a[i*n + j] * x[j];
		//std::cout << sum << " " << b[i] << std::endl;
		if (abs(sum - a[i*n + n + 1]) >= EPS) {
			verdict = false;
			//std::cout << i << ". " << sum << " " << b[i] << " " << abs(sum - b[i]) << std::endl;
		}
	}

	return verdict;
}

int main()
{
	int k = 0, p = 2, i = 0, j = 0, d = 0, d_i = 0;
	int n, thread_count;
	cin >> n >> thread_count;
	omp_set_dynamic(false);
	omp_set_num_threads(thread_count);

	double *a = (double*)malloc(sizeof(double) * n * (n + 1));
	double *ans = (double*)malloc(sizeof(double) * n);

	init(a, ans, n);

	double finish_time, start_time;

	start_time = omp_get_wtime();
#pragma omp parallel for private(k, i, j, d, p) shared(a) num_threads(n_threads)
	for (k = 0; k < n; ++k) {
		p = 2;
		i = k;
		d = get_deep(n - k - 1); //(n - k - 1) % 2==0 ? n - k - 1 : n - k;
		d_i = 0;
		for (; n - k >= p / 2; i += p, ++d_i) {
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

				a[k*n + j]				 = c * old_value_a_kj + s * old_value_a_ij;
				a[(i - k / 2 - 1)*n + j] = c * old_value_a_ij - s * old_value_a_kj;
			}

			if (i >= n - 2)
				p *= 2, i = k;
		}
	}
	finish_time = omp_get_wtime();

	int b_index = n + 1;
	long double sum = 0;
	i = n, j = n;
#pragma omp parallel private (i, j) reduction(+:sum) num_threads(*thread_count/*omp_get_max_threads()*/)
	{
#pragma omp for schedule(dynamic)
		for (i = n; i > 0; --i) {
			sum = 0;
			for (j = n; j > i; --j)
				sum += ans[j] * a[i*n + j];

			ans[i] = (a[i*n + b_index] - sum) / a[i*n + i];
		}
	}

	printf("Ax=b? ");
	if (check(a, ans, n) == 1)
		printf("Yes\n");
	else
		printf("No\n");
	printf("TimeOmp = %f\n", double(finish_time - start_time));
	//print_result(a, n);

	//getch();
	return 0;
}


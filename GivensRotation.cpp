/*
Steshenko A.S.
*/

#include <stdio.h>
#include <conio.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <set>
#include <map>
#include <stack>
#include <queue>
#include <deque>
#include <string>
#include <list>
#include <ctime>
#include <cctype>

#define rep(i, j, n) for(int i = j;i <= n; i++)
#define rev(i, j, n) for(int i = n;i >= j; i--)

//using namespace std;  // don`t use namespace (it`s bad)

int main()
{
	time_t START = clock();
	std::ios_base::sync_with_stdio(0);
	std::cin.tie(0);
	
	int n;
	std::cin >> n;

	std::vector <double> x(n + 1);
	std::vector <std::vector <double> > a(n + 1);

	for (int i = 0; i <= n; ++i) a[i].assign(n + 2, 0);

	rep(i, 1, n)
		rep(j, 1, n + 1)
			std::cin >> a[i][j];

	/*
	rep(i, 1, n) {
		rep(j, 1, n + 1)
			std::cout << a[i][j] << " ";
		std::cout << endl;
	}
	*/

	int b_index = n + 1; // matrix A and vector b is lie in single matrix A
	// start clock
	for (int k = 1; k < n; ++k) {
		for (int i = k + 1; i <= n; ++i) {
			double  den = pow(a[k][k] * a[k][k] + a[i][k] * a[i][k], 0.5);
			double  c = a[k][k] / den, 
					s = a[i][k] / den;

			for (int j = k; j <= n; ++j) {
				double  old_value_a_ij = a[i][j],
						old_value_a_kj = a[k][j];

				a[k][j] = c * old_value_a_kj + s * old_value_a_ij;
				a[i][j] = c * old_value_a_ij - s * old_value_a_kj;
			}

			double  old_value_a_kb = a[k][b_index],
					old_value_a_ib = a[i][b_index];

			a[k][b_index] = c * old_value_a_kb + s * old_value_a_ib;
			a[i][b_index] = c * old_value_a_ib - s * old_value_a_kb;
		}
	}

	x[n] = a[n][b_index] / a[n][n];
	for (int i = n - 1; i > 0; --i) {
		double sum = 0;
		for (int j = n; j > i; --j) sum += x[j] * a[i][j];

		x[i] = (a[i][b_index] - sum) / a[i][i];
	}
	// finish clock

	/*
	rep(i, 1, n) {
		rep(j, 1, n + 1)
			std::cout << a[i][j] << " ";
		std::cout << endl;
	}
	*/

	for (int i = 1; i < x.size(); ++i)
		std::cout << "x[" << i << "] = " << x[i] << std::endl;

	time_t FINISH = clock();
	std::cerr << "Time = " << double(FINISH - START) / CLOCKS_PER_SEC << std::endl;

	getch();
	return 0;
}


#include <cassert>
#include <climits>
#include <cmath>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <numeric>
#include <thread>
#include <vector>

#include <conversions.h>
#include <matrix_helper.h>
#include <matrix_operations.h>
#include <signals.h>
#include <vector_templates.h>

#include "enumeration.h"

using namespace NTL;
using namespace std;

static unsigned int number_of_max_threads = thread::hardware_concurrency();
static unsigned int number_of_threads = 0;
static mutex mut_thread_count;
static mutex mut_cout;
static condition_variable is_thread_free;

// could be replaced by atomic<double>
static mutex mut_errL2;
static double current_errL2 = numeric_limits<double>::max();

double LengthPruningOptThread(shared_ptr<vector<long>> p_solution,
							  matrix<long> const &B, matrix<double> const &mu,
							  vector<double> const &R, vector<double> const &t,
							  vector<long> c, size_t lvl,
							  vector<double> errVec);

vector<long> LengthPruningOpt(matrix<long> const &B,
							  vector<double> const &rbound,
							  vector<double> const &t, long q,
							  matrix<double> const &mu) {
	// initialize
	size_t m = B.size();
	matrix<double> sigma(m + 1, vector<double>(m, 0));
	vector<double> r(m + 1);
	for (size_t i = 0; i < m + 1; ++i)
		r[i] = i;
	vector<double> c(m);
	vector<long> v(m);
	vector<long> w(m);
	vector<long> v1(m);
	vector<double> rho(m + 1, 0);

	double current_errL = 0.0;

	vector<long> sol(m, 0); // solution-vector --EK

	// babai
	for (size_t k = m; k > 0; --k) {
		for (size_t i = m; i > k; --i) {
			sigma[i - 1][k - 1] =
				sigma[i][k - 1] + (t[i - 1] - v[i - 1]) * mu[i - 1][k - 1];
		}
		c[k - 1] = t[k - 1] + sigma[k][k - 1];
		v[k - 1] = lround(c[k - 1]);
		w[k - 1] = 1;
		rho[k - 1] = rho[k] + pow(c[k - 1] - v[k - 1], 2) * mu[k - 1][k - 1];
	}

	v1 = v;
	size_t k = 0;
	current_errL = rho[k];

	while (true) {
		rho[k] = rho[k + 1] + pow((c[k] - v[k]), 2) * mu[k][k];
		if (rho[k] <= rbound[k]) // we generate the sequence of the form R_1^2 >
								 // ... > R_m^2
		{
			if (k == 0) {
				if (rho[k] <= current_errL) {
					cout << "new solution found " << v << endl;
					// save v:
					v1 = v;
					current_errL = rho[k];
					cout << "of length " << rho[k] << endl;
				}

				goto TRAVERSE_UP;
			} else {
				k--;
				r[k] = max(r[k], r[k + 1]);
				for (size_t i = r[k]; i > k; --i)
					sigma[i][k] = sigma[i + 1][k] + (t[i] - v[i]) * mu[i][k];
				c[k] = t[k] + sigma[k + 1][k];
				v[k] = lround(c[k]);
				w[k] = 1;
			}
		} else {
		TRAVERSE_UP:
			k++;
			if (k == m) {
				cout << "no new solution found " << v1 << endl;
				for (size_t i = 0; i < m; ++i) {
					sol += v1[i] * B[i];
				}
				cout << "sol = " << sol % q << endl;
				return sol % q;
			}

			r[k] = k;
			if (v[k] > c[k])
				v[k] -= w[k];
			else
				v[k] += w[k];
			w[k]++;
		}
	}
}

vector<long> LengthPruningOptParall(matrix<long> const &B,
									matrix<double> const &B_star,
									vector<double> const &rbound,
									vector<long> const &t, long q, size_t lvl) {
	size_t spawned_threads = 0;
	// initialize
	const matrix<double> mu = muGSO(B);
	vector<double> mu_t = muT(B_star, t);

	size_t j = 0;
	size_t m = B.size();
	vector<double> t_coeff = coeffs(B, t);
	vector<double> c_star(lvl);
	vector<long> c_min(lvl);
	vector<long> c_max(lvl);
	vector<long> c(m);
	vector<future<double>> future_solutions;
	vector<shared_ptr<vector<long>>> solutions;

	vector<long> t_err(t);
	vector<long> solution(m);
	vector<double> errVec(m + 1, 0);

	matrix<double> sigma(m + 1, vector<double>(m, 0));
	sigma[m] = mu_t;

	// iterate down to this level and spawn threads
	while (!got_sigterm) {
		// update c_min[j], c_max[j]
		c_star[j] = sigma[m - j][m - j - 1];
		double interval = sqrt(fabs(rbound[m - j - 1] - errVec[j]) /
							   mu[m - j - 1][m - j - 1]);
		c_min[j] = (long)ceil(c_star[j] - interval / 2.0);
		c_max[j] = (long)floor(c_star[j] + interval / 2.0);
		assert((lround(c_star[j]) >= c_min[j]) &&
			   (lround(c_star[j]) <= c_max[j]) &&
			   "c_star not in the interval!");

		// TRAVERSE DOWN
		c[j] = c_min[j];

		if (j != lvl - 1) {
			errVec[j + 1] =
				errVec[j] + pow(c_star[j] - c[j], 2) * mu[m - j - 1][m - j - 1];

			for (long i = m - j - 1; i >= 0; --i)
				sigma[m - j - 1][i] =
					sigma[m - j][i] - (double)c[j] * mu[m - j - 1][i];

			t_coeff[m - j - 1] -= c[j];

			j++;
		} else {
			// main loop to start threads
			while (c[j] <= c_max[j]) {
				vector<double> errVec_new = errVec;
				errVec_new[j + 1] =
					errVec[j] +
					pow(c_star[j] - c[j], 2) * mu[m - j - 1][m - j - 1];

				t_coeff[m - j - 1] -= c[j];

				shared_ptr<vector<long>> p_solution(new vector<long>(m, 0));
				for (size_t i = 0; i <= j; ++i) {
					*p_solution += c[i] * B[m - i - 1];
					*p_solution %= q;
				}
				solutions.push_back(p_solution);

				{ // lock mutex for thread count, wait until we can spawn
					// another thread and spawn new thread
					unique_lock<mutex> lock{mut_thread_count};
					is_thread_free.wait(lock, [&] {
						return number_of_threads < number_of_max_threads;
					});
					// we tell the thread, to start at level lvl+1 here,
					// because we already have #lvl (0 to lvl-1) iterations
					// done before and are doing the lvl+1'th iteration in
					// the current while loop
					future_solutions.push_back(async(
						launch::async, LengthPruningOptThread, solutions.back(),
						B, mu, rbound, t_coeff, c, j + 1, errVec_new));
					number_of_threads++;
				} // mutex is freed, when lock goes out of scope
				is_thread_free.notify_one();
				spawned_threads++;
#ifdef DEBUG
				{
					unique_lock<mutex> lock{mut_cout};
					cout << "spawned thread no. " << spawned_threads << endl;
				}
#endif

				c[j]++;
			}

			// 2. traverse up until there is at least one sibbling to the right
			do {
				if (j == 0 && c[j] >= c_max[j])
					goto RETURN_SPAWN_THREADS;
				t_coeff[m - j - 1] += c[j];
				j--;
			} while (c[j] >= c_max[j]);

			// TRAVERSE RIGHT
			c[j]++;
			errVec[j + 1] +=
				(-2.0 * c_star[j] + 2.0 * c[j] - 1) * mu[m - j - 1][m - j - 1];

			for (long i = m - j - 1; i >= 0; --i)
				sigma[m - j - 1][i] -= mu[m - j - 1][i];

			t_coeff[m - j - 1]--;

			j++;
		}
	}

RETURN_SPAWN_THREADS:
	// after every thread is started, gather the results
	auto actual_solution = solutions.begin();
	vector<long> &current_solution = **actual_solution;
	for (auto &future_errL2 : future_solutions) {
		// get blocks until the thread has finished
		auto const &actual_errL2 = future_errL2.get();

		// 1. check wether the leaf corresponds to a new candidate solution
		//    (distance of actual solution to original input equals length
		//     of actual t vector, while solution is coded by the c vector)
		if (actual_errL2 <= current_errL2) {
			{ // lock mutex for errL2
				unique_lock<mutex> lock{mut_errL2};
				current_errL2 = actual_errL2;
			} // mutex is freed
			current_solution = **actual_solution;
			// assert(current_solution.size() == t.size() &&
			//		"solution has different size than t");
		}
		actual_solution++;
	}

	{
		unique_lock<mutex> lock{mut_cout};
		cout << "overall spawned threads: " << spawned_threads << endl;
	}
	// return the solution from the thread
	return current_solution % q;
}

double LengthPruningOptThread(shared_ptr<vector<long>> p_solution,
							  matrix<long> const &B, matrix<double> const &mu,
							  vector<double> const &rbound,
							  vector<double> const &t, vector<long> v,
							  size_t lvl, vector<double> rho) {
	// initialize
	// size_t m = B.size()-lvl;
	size_t m = B.size();
	vector<double> r(m + 1);
	for (size_t i = 0; i < m + 1; ++i)
		r[i] = i;
	vector<double> c(m);
	vector<long> w(m);
	vector<long> v1(m);
	// vector<double> rho(m+1, 0);

	matrix<double> sigma(m + 1, vector<double>(m, 0));

	double thread_errL2 = 0.0;

	vector<long> sol(m, 0); // solution-vector --EK

	// babai
	for (size_t k = m; k > 0; --k) {
		for (size_t i = m; i > k; --i)
			sigma[i - 1][k - 1] =
				sigma[i][k - 1] + (t[i - 1] - v[i - 1]) * mu[i - 1][k - 1];
		c[k - 1] = t[k - 1] + sigma[k][k - 1];
		v[k - 1] = lround(c[k - 1]);
		w[k - 1] = 1;
		rho[k - 1] = rho[k] + pow(c[k - 1] - v[k - 1], 2) * mu[k - 1][k - 1];
	}

	v1 = v;
	size_t k = 0;
	thread_errL2 = rho[k];

	m = m - lvl;

	while (true) {
		rho[k] = rho[k + 1] + pow((c[k] - v[k]), 2) * mu[k][k];
		if (rho[k] <= rbound[k]) // we generate the sequence of the form R_1^2 >
								 // ... > R_m^2
		{
			if (k == 0) {
				if (rho[k] <= thread_errL2) {
					cout << "new solution found " << v << endl;
					// save v:
					v1 = v;
					thread_errL2 = rho[k];
					cout << "of length " << rho[k] << endl;
				}

				goto TRAVERSE_UP;
			} else {
				k--;
				r[k] = max(r[k], r[k + 1]);
				for (size_t i = r[k]; i > k; --i)
					sigma[i][k] = sigma[i + 1][k] + (t[i] - v[i]) * mu[i][k];
				c[k] = t[k] + sigma[k + 1][k];
				v[k] = lround(c[k]);
				w[k] = 1;
			}
		} else {
		TRAVERSE_UP:
			k++;
			if (k == m) {
				for (size_t i = 0; i < m; ++i)
					*p_solution += v1[i] * B[i];
				cout << "sol = " << *p_solution << endl;
				goto RETURN_LENGTH_THREAD;
			}

			r[k] = k;
			if (v[k] > c[k])
				v[k] -= w[k];
			else
				v[k] += w[k];
			w[k]++;
		}
	}

RETURN_LENGTH_THREAD : { // lock mutex for thread count
	unique_lock<mutex> lock{mut_thread_count};
	number_of_threads--;
} // mutex is freed, when lock goes out of scope
	is_thread_free.notify_one();

	return thread_errL2;
}

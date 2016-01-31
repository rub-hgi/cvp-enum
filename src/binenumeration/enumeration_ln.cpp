#include <atomic>
#include <cassert>
#include <climits>
#include <cmath>
#include <condition_variable>
#include <deque>
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
static mutex mut_err;
static double current_error = numeric_limits<double>::max();

tuple<vector<long>, double>
LengthPruningOptThread(matrix<long> const &B, matrix<double> const &mu,
					   vector<double> const &b_star_lengths,
					   vector<double> const &rbound, vector<double> const &mu_t,
					   size_t lvl, double rho_init, long q,
					   vector<long> v_breadth);

vector<long> LengthPruningOpt(matrix<long> const &B,
							  matrix<double> const &B_star,
							  vector<double> const &rbound,
							  vector<long> const &t, long q) {
	// initialize
	const matrix<double> mu = muGSO(B);
	const vector<double> b_star_lengths = VectorLengths(B_star);
	vector<double> mu_t = muT(B_star, t);

	size_t m = B.size();
	matrix<double> sigma(m + 1, vector<double>(m, 0));
	sigma[m] = mu_t;
	vector<double> c(m);
	vector<long> w(m);
	vector<long> v(m);

	vector<double> rho(m + 1, 0);

	// babai for mu_t
	for (size_t k = m; k > 0; --k) {
		c[k - 1] = sigma[k][k - 1];
		v[k - 1] = lround(c[k - 1]);
		w[k - 1] = 1;
		rho[k - 1] =
			rho[k] + pow(c[k - 1] - v[k - 1], 2) * b_star_lengths[k - 1];
		for (size_t i = 0; i < k; i++)
			sigma[k - 1][i] = sigma[k][i] - (double)v[k - 1] * mu[k - 1][i];
	}

	vector<long> v1(v);
	size_t k = 0;
	double error = rho[k];

	cout << "finished with babai's part, start while loop" << endl;
	while (!got_sigterm) {
		rho[k] = rho[k + 1] + pow((c[k] - v[k]), 2) * b_star_lengths[k];
		if (rho[k] <= rbound[k]) // we generate the sequence of the form R_1^2 >
								 // ... > R_m^2
		{
			if (k == 0) {
				if (rho[k] <= error) {
					v1 = v;
					error = rho[k];
#ifdef DEBUG
					cout << "new solution of length " << rho[k] << " found"
						 << endl;
					cout << "v         = " << v << endl;
#endif
				}

				goto TRAVERSE_UP_LENGTH;
			} else {
				for (size_t i = 0; i < k + 1; ++i)
					sigma[k][i] = sigma[k + 1][i] - double(v[k]) * mu[k][i];
				k--;
				c[k] = sigma[k + 1][k];
				v[k] = lround(c[k]);
				w[k] = 1;
			}
		} else {
		TRAVERSE_UP_LENGTH:
			k++;
			if (k == m)
				break;

			if (v[k] > c[k])
				v[k] -= w[k];
			else
				v[k] += w[k];
			w[k]++;
		}
	}

	vector<long> sol(B.size(), 0);
	for (size_t i = 0; i < m; ++i) {
		sol += v1[i] * B[i];
	}
#ifdef DEBUG
	cout << "v   = " << v1 << endl;
	cout << "sol = " << sol % q << endl;
#endif
	return sol % q;
}

typedef tuple<vector<double>, // mu_t
			  double,		  // err
			  vector<long>	// v
			  > thread_args;

vector<long> LengthPruningOptParall(matrix<long> const &B,
									matrix<double> const &B_star,
									vector<double> const &rbound,
									vector<long> const &t, long q,
									long n_threads, long factor_lvl) {
	if (n_threads != 0) {
		number_of_max_threads = n_threads;
	} else {
		n_threads = number_of_max_threads;
	}

	// step 1: Algorithm 3
	// Traverse Breadth-First

	// initialize
	const size_t m = B.size();
	deque<thread_args> queue;
	vector<long> v(m, 0);
	queue.push_back(thread_args(vector<double>(muT(B_star, t)), 0.0, v));

	const matrix<double> mu = muGSO(B);
	const vector<double> b_star_lengths = VectorLengths(B_star);

	long n = 1;
	long k = m - 1;

	// factor_lvl balances short running threads
	while ((n < factor_lvl * n_threads) && k > -1) {
		long cnt = n;
		n = 0;
		for (long i = 1; i <= cnt; ++i) {
			auto nodeargs = queue.front();
			// TODO optimization: move vector instead of copy it?
			vector<double> mu_t = get<0>(nodeargs);
			double error = get<1>(nodeargs);
			vector<long> v = get<2>(nodeargs);
			queue.pop_front();

			long interval =
				ceil(sqrt(fabs(rbound[k] - error) / b_star_lengths[k]));
			// n += interval;
			double c_star = mu_t[k];
			for (long j = 0; j < interval; ++j) {
				v[k] = lround(c_star + j);
				queue.push_back(thread_args(
					vector<double>(mu_t - (double)v[k] * mu[k]),
					error + pow((c_star - v[k]), 2) * b_star_lengths[k], v));
				n++;

				assert(!((lround(c_star - j) == v[k]) && (j != 0)) &&
					   "whoops, this wasn't expected");
				if (lround(c_star - j) != v[k]) {
					v[k] = lround(c_star - j);
					queue.push_back(thread_args(
						vector<double>(mu_t - (double)v[k] * mu[k]),
						error + pow((c_star - v[k]), 2) * b_star_lengths[k],
						v));
					n++;
				}
			}

			// corner-cases
			if (error +
					pow(c_star - lround(c_star) + interval, 2) *
						b_star_lengths[k] <
				rbound[k]) {
				v[k] = lround(c_star) - interval;
				queue.push_back(thread_args(
					vector<double>(mu_t - (double)v[k] * mu[k]),
					error + pow((c_star - v[k]), 2) * b_star_lengths[k], v));
				n++;
			}

			if (error +
					pow(c_star - lround(c_star) - interval, 2) *
						b_star_lengths[k] <
				rbound[k]) {
				v[k] = lround(c_star) + interval;
				queue.push_back(thread_args(
					vector<double>(mu_t - (double)v[k] * mu[k]),
					error + pow((c_star - v[k]), 2) * b_star_lengths[k], v));
				n++;
			}
		}
		k--;
	}
	cout << "breadth-first traversal till k = " << k << endl;

	// check if Babai's only is enough
	if (k < 0) {
		cout << "finished with babai, traversed whole tree, computing solution"
			 << endl;
		vector<double> mu_t;
		double error = numeric_limits<double>::max();
		while (!queue.empty()) {
			auto nodeargs = queue.front();
			queue.pop_front();
			if (error > get<1>(nodeargs)) {
				// TODO optimization: move vector instead of copy it?
				mu_t = get<0>(nodeargs);
				error = get<1>(nodeargs);
			}
		}
		vector<double> err_prime(m, 0);
		for (size_t i = 0; i < m; ++i)
			err_prime += mu_t[i] * B_star[i];
		round(err_prime);
		return (t - to_stl<long>(err_prime)) % q;
	}

// step 2
// call LengthPruningThreads for every item in queue
#ifdef DEBUG
	size_t spawned_threads = 0;
#endif
	vector<future<tuple<vector<long>, double>>> future_solutions;
	while (!queue.empty()) {
		auto nodeargs = queue.front();
		queue.pop_front();
		// TODO optimization: move vector instead of copy it?
		vector<double> mu_t = get<0>(nodeargs);
		double error = get<1>(nodeargs);
		vector<long> v = get<2>(nodeargs);
		{ // lock mutex for thread count, wait until we can spawn
			// another thread and spawn new thread
			unique_lock<mutex> lock{mut_thread_count};
			is_thread_free.wait(lock, [&] {
				return number_of_threads < number_of_max_threads;
			});

			future_solutions.push_back(
				async(launch::async, LengthPruningOptThread, B, mu,
					  b_star_lengths, rbound, mu_t, k, error, q, v));

			number_of_threads++;
		} // mutex is freed, when lock goes out of scope
#ifdef DEBUG
		spawned_threads++;
		{
			unique_lock<mutex> lock{mut_cout};
			cout << "spawned thread no. " << spawned_threads << endl;
		}
#endif
	}

	// step 3
	// collect best solution from threads
	vector<long> current_solution(m, 0);
	for (auto &future_solution : future_solutions) {
		auto const &solution = future_solution.get();
		// TODO optimization: move vector instead of copy it?
		vector<long> sol = get<0>(solution);
		double error = get<1>(solution);

		if (error <= current_error) {
			{ // lock mutex for current_error
				unique_lock<mutex> lock{mut_err};
				current_error = error;
			}
			// TODO optimization: move vector instead of copy it?
			current_solution = sol;
		}
	}
	return current_solution % q;
}

tuple<vector<long>, double>
LengthPruningOptThread(matrix<long> const &B, matrix<double> const &mu,
					   vector<double> const &b_star_lengths,
					   vector<double> const &rbound, vector<double> const &mu_t,
					   size_t lvl, double rho_init, long q,
					   vector<long> v_breadth) {
	// initialize
	size_t m = lvl + 1; // TODO correct?
	matrix<double> sigma(m + 1, vector<double>(B.size(), 0));
	sigma[m] = mu_t;
	vector<double> c(B.size());
	vector<long> w(m);
	vector<long> v(m);

	vector<double> rho(m + 1, 0);
	rho[m] = rho_init;

	// babai for mu_t
	for (size_t k = m; k > 0; --k) {
		c[k - 1] = sigma[k][k - 1];
		v[k - 1] = lround(c[k - 1]);
		w[k - 1] = 1;
		rho[k - 1] =
			rho[k] + pow(c[k - 1] - v[k - 1], 2) * b_star_lengths[k - 1];
		for (size_t i = 0; i < k; i++)
			sigma[k - 1][i] = sigma[k][i] - (double)v[k - 1] * mu[k - 1][i];
	}

	vector<long> v1(v);
	size_t k = 0;
	double error = rho[k];

	while (!got_sigterm) {
		rho[k] = rho[k + 1] + pow((c[k] - v[k]), 2) * b_star_lengths[k];
		if (rho[k] <= rbound[k]) // we generate the sequence of the form R_1^2 >
								 // ... > R_m^2
		{
			if (k == 0) {
				if (rho[k] <= current_error) {
					v1 = v;
					error = rho[k];
					{ // lock mutex for errL2
						unique_lock<mutex> lock{mut_err};
						current_error = rho[k];
					}
#ifdef DEBUG
					{ // lock mutex for cout
						unique_lock<mutex> lock{mut_cout};
						cout << "new solution of length " << rho[k] << " found"
							 << endl;
						cout << "v         = " << v << endl;
						cout << "v_breadth = " << v_breadth << endl;
					}
#endif
				}

				goto TRAVERSE_UP_LENGTH_THREAD;
			} else {
				for (size_t i = 0; i < k + 1; ++i)
					sigma[k][i] = sigma[k + 1][i] - (double)v[k] * mu[k][i];
				k--;
				c[k] = sigma[k + 1][k];
				v[k] = lround(c[k]);
				w[k] = 1;
			}
		} else {
		TRAVERSE_UP_LENGTH_THREAD:
			k++;
			if (k == m)
				break;

			if (v[k] > c[k])
				v[k] -= w[k];
			else
				v[k] += w[k];
			w[k]++;
		}
	}

	vector<long> sol(B.size(), 0);
	for (size_t i = 0; i < m; ++i)
		sol += v1[i] * B[i];
	for (size_t i = m; i < B.size(); i++)
		sol += v_breadth[i] * B[i];

	{ // lock mutex for thread count
		unique_lock<mutex> lock{mut_thread_count};
		number_of_threads--;
	}
	is_thread_free.notify_one();

	return make_tuple(sol, error);
}

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
#include <matrix_operations.h>
#include <signals.h>
#include <vector_templates.h>

#include "enumeration.h"

using namespace std;

static unsigned int number_of_max_threads = thread::hardware_concurrency();
static unsigned int number_of_threads = 0;
static mutex mut_thread_count;
static mutex mut_cout;
static condition_variable is_thread_free;

// could be replaced by atomic<double>
static mutex mut_errL2;
static double current_errL2 = numeric_limits<double>::max();

double NearestPlanesLPOptThread(
	shared_ptr<vector<long>> p_solution,
	matrix<long> const& B,
	matrix<double> const& mu,
	vector<long> const& d,
	matrix<double> sigma,
	vector<long> c,
	size_t j_new,
	vector<double> errVec
);

vector<long> NearestPlanesLPOpt(matrix<long> const& B, vector<long> const& d,
		vector<long> const& t, long q, matrix<double> const& B_star)
{
	// initialize
	size_t j = 0;
	size_t m = B.size();
	vector<double> c_star(m);
	vector<long> c_min(m);
	vector<long> c_max(m);
	vector<long> c(m);

	matrix<double> mu = muGSO(B);
	vector<double> mu_t = muT(B_star, t);
	vector<long> t_err(t);
	vector<long> solution(m);
	vector<double> errVec(m+1, 0);

	matrix<double> sigma(m+1, vector<double>(m, 0));
	double current_min = numeric_limits<double>::max();
	sigma[m] = mu_t;

	// traverse search tree depth first
	while (!got_sigterm)
	{
		if (j != m)
		{
			// update c_min[j], c_max[j]
			c_star[j] = sigma[m-j][m-j-1];
			c_min[j] = ceil(c_star[j] - d[m-j-1] / 2.0);
			c_max[j] = floor(c_star[j] + d[m-j-1] / 2.0);
			//assert( abs(c_max[j] - c_min[j]) == d[m-j-1] - 1
			//		&& "c_min/max intervall too big" );
			assert((lround(c_star[j]) >= c_min[j])
				&& (lround(c_star[j]) <= c_max[j])
				&& "c_star not in the interval!" );

			// TRAVERSE DOWN
			c[j] = c_min[j];

			errVec[j+1] = errVec[j] + pow(c_star[j]-c[j], 2) * mu[m-j-1][m-j-1];

			for (long i=m-j-1; i>=0; --i)
				sigma[m-j-1][i] = sigma[m-j][i] - (double) c[j] * mu[m-j-1][i];

			j++;
			if (errVec[j] > current_min)
				goto TRAVERSE_UP_LP;

		} else
		{
			if (errVec[j] < current_min)
			{
				current_min = errVec[j];
				solution = c;
			}

TRAVERSE_UP_LP:
			// 2. traverse up until there is at least one sibbling to the right
			do {
				if (j == 0 && c[j] >= c_max[j])
					goto RETURN_LP;
				j--;
			} while (c[j] >= c_max[j]);
			// TRAVERSE RIGHT
			c[j]++;

			errVec[j+1] += (-2.0 * c_star[j] + 2.0 * c[j] - 1) * mu[m-j-1][m-j-1];

			for (long i=m-j-1; i >= 0; --i)
				sigma[m-j-1][i] -= mu[m-j-1][i];

			j++;
		}
	}

RETURN_LP:
	vector<long> sol(m, 0);
	cout << " coefficients: " << solution << endl;
	for (size_t i=0; i<B.size(); ++i)
		sol += solution[i] * B[m-i-1];
	return sol % q;
}

vector<long> NearestPlanesLPOptParall(matrix<long> const& B,
		matrix<double> const& B_star, vector<long> const& d, vector<long> t, long q,
		size_t lvl)
{
	size_t spawned_threads = 0;
	// initialize
	const matrix<double> mu = muGSO(B);
	vector<double> mu_t = muT(B_star, t);

	size_t j = 0;
	size_t m = B.size();
	vector<double> c_star(lvl);
	vector<long> c_min(lvl);
	vector<long> c_max(lvl);
	vector<long> c(m);
	vector<future<double>> future_solutions;
	vector<shared_ptr<vector<long>>> solutions;

	vector<long> t_err(t);
	vector<long> solution(m);
	vector<double> errVec(m+1, 0);

	matrix<double> sigma(m+1, vector<double>(m, 0));
	sigma[m] = mu_t;

	// iterate down to this level and spawn threads
	while (!got_sigterm)
	{
		// update c_min[j], c_max[j]
		c_star[j] = sigma[m-j][m-j-1];
		c_min[j] = (long) ceil(c_star[j] - (double) d[m-j-1] / 2.0);
		c_max[j] = (long) floor(c_star[j] + (double) d[m-j-1] / 2.0);
		assert(abs(c_max[j] - c_min[j]) == d[m-j-1] - 1
			&& "c_min/max intervall too big" );
		assert((lround(c_star[j]) >= c_min[j])
			&& (lround(c_star[j]) <= c_max[j])
			&& "c_star not in the interval!" );

		// TRAVERSE DOWN
		c[j] = c_min[j];

		if (j != lvl-1)
		{
			errVec[j+1] = errVec[j] + pow(c_star[j]-c[j], 2) * mu[m-j-1][m-j-1];

			for (long i=m-j-1; i>=0; --i)
				sigma[m-j-1][i] = sigma[m-j][i] - (double) c[j] * mu[m-j-1][i];

			j++;
		} else
		{
			// main loop to start threads
			while (c[j] <= c_max[j])
			{
				vector<double> errVec_new = errVec;
				errVec_new[j+1] = errVec[j] + pow(c_star[j]-c[j], 2) * mu[m-j-1][m-j-1];
				matrix<double> sigma_new = sigma;
				// current solution for new thread
				for (long i=m-j-1; i>=0; --i)
					sigma_new[m-j-1][i] = sigma[m-j][i] - (double) c[j] * mu[m-j-1][i];

				shared_ptr<vector<long>> p_solution(new vector<long>(m, 0));
				for (size_t i = 0; i <= j; ++i)
				{
					*p_solution += c[i] * B[m-i-1];
					*p_solution %= q;
				}
				solutions.push_back(p_solution);

				{   // lock mutex for thread count, wait until we can spawn
					// another thread and spawn new thread
					unique_lock<mutex> lock { mut_thread_count };
					is_thread_free.wait(lock, [&]{
						return number_of_threads < number_of_max_threads;
						});
					// we tell the thread, to start at level lvl+1 here,
					// because we already have #lvl (0 to lvl-1) iterations
					// done before and are doing the lvl+1'th iteration in
					// the current while loop
					future_solutions.push_back(async(launch::async,
						NearestPlanesLPOptThread, solutions.back(), B, mu, d, sigma_new,
						c, j+1, errVec_new));
					number_of_threads++;
				}	// mutex is freed, when lock goes out of scope
				is_thread_free.notify_one();
				spawned_threads++;

				c[j]++;
			}

			// 2. traverse up until there is at least one sibbling to the right
			do {
				if (j == 0 && c[j] >= c_max[j])
					goto RETURN_SPAWN_THREADS;
				j--;
			} while (c[j] >= c_max[j]);

			// TRAVERSE RIGHT
			c[j]++;
			errVec[j+1] += (-2.0 * c_star[j] + 2.0 * c[j] - 1) * mu[m-j-1][m-j-1];

			for (long i=m-j-1; i >= 0; --i)
				sigma[m-j-1][i] -= mu[m-j-1][i];

			j++;
		}
	}

RETURN_SPAWN_THREADS:
	// after every thread is started, gather the results
	auto actual_solution = solutions.begin();
	vector<long> &current_solution = **actual_solution;
	for (auto &future_errL2 : future_solutions)
	{
		// get blocks until the thread has finished
		auto const& actual_errL2 = future_errL2.get();

		// 1. check wether the leaf corresponds to a new candidate solution
		//    (distance of actual solution to original input equals length
		//     of actual t vector, while solution is coded by the c vector)
		if (actual_errL2 <= current_errL2)
		{
			{   // lock mutex for errL2
				unique_lock<mutex> lock { mut_errL2 };
				current_errL2 = actual_errL2;
			}	// mutex is freed
			current_solution = **actual_solution;
			//assert(current_solution.size() == t.size() &&
			//		"solution has different size than t");
		}
		actual_solution++;
	}

	{
		unique_lock<mutex> lock { mut_cout };
		cout << "overall spawned threads: " << spawned_threads << endl;
	}
	// compute the solution from the returned c vector
	vector<long> sol(m, 0);
	for (size_t i = 0; i < m; ++i)
		sol += current_solution[i] * B[m-i-1];
	return sol % q;
}

double NearestPlanesLPOptThread(
	shared_ptr<vector<long>> p_solution,
	matrix<long> const& B,
	matrix<double> const& mu,
	vector<long> const& d,
	matrix<double> sigma,
	vector<long> c,
	size_t lvl,
	vector<double> errVec)
{
	// initialize
	size_t j = lvl;
	size_t m = B.size();
	vector<double> c_star(m);
	vector<long> c_min(m);
	vector<long> c_max(m);

	double thread_min = numeric_limits<double>::max();

	// traverse search tree depth first
	while (!got_sigterm)
	{
		if (j != m)
		{
			// update c_min[j], c_max[j]
			c_star[j] = sigma[m-j][m-j-1];
			c_min[j] = ceil(c_star[j] - d[m-j-1] / 2.0);
			c_max[j] = floor(c_star[j] + d[m-j-1] / 2.0);
			//assert( abs(c_max[j] - c_min[j]) == d[m-j-1] - 1
			//		&& "c_min/max intervall too big" );
			assert((lround(c_star[j]) >= c_min[j])
					&& (lround(c_star[j]) <= c_max[j])
					&& "c_star not in the interval!" );

			// TRAVERSE DOWN
			c[j] = c_min[j];
			errVec[j+1] = errVec[j] + pow(c_star[j]-c[j], 2) * mu[m-j-1][m-j-1];

			for (long i=m-j-1; i>=0; --i)
				sigma[m-j-1][i] = sigma[m-j][i] - (double) c[j] * mu[m-j-1][i];

			j++;
			if (errVec[j] > current_errL2)
				goto TRAVERSE_UP_LP_THREAD;
		} else
		{
			// 1. check wether the leaf corresponds to a new candidate solution
			//    (distance of actual solution to original input equals length
			//     of actual t vector, while solution is coded by the c vector)
			if (errVec[j] < current_errL2)
			{
				{   // lock mutex for errL2
					unique_lock<mutex> lock { mut_errL2 };
					current_errL2 = errVec[j];
				}	// mutex is freed
				thread_min = errVec[j];
				// return c vector as solution, that means, we have to compute
				// the solution in the main function, but thats ok - and should even
				// save some time, as we only need to compute the final solution
				*p_solution = c;
			}
TRAVERSE_UP_LP_THREAD:
			// 2. traverse up until there is at least one sibbling to the right
			do {
				assert(j >= lvl && "j < lvl, traversed too high");
				if (j == lvl && c[j] >= c_max[j])
					// if no siblings left, break the loop,
					// decrement the number of active threads and
					// return the current solution's error
					goto RETURN_LP_THREAD;
				j--;
			} while (c[j] >= c_max[j]);

			// TRAVERSE RIGHT
			c[j]++;
			errVec[j+1] += (-2.0 * c_star[j] + 2.0 * c[j] - 1) * mu[m-j-1][m-j-1];

			for (long i=m-j-1; i >= 0; --i)
				sigma[m-j-1][i] -= mu[m-j-1][i];

			j++;
		}
	}
RETURN_LP_THREAD:
	{   // lock mutex for thread count
		unique_lock<mutex> lock { mut_thread_count };
		number_of_threads--;
	}	// mutex is freed, when lock goes out of scope
	is_thread_free.notify_one();

	return thread_min;
}


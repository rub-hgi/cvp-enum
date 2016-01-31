// ----------------------------------------------------------------------------
// Title      : Vector Templates
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : vector_templates.h
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Templated implementations for vector operations. This needs to be
//     in a header file, as c++ does not support templated implementations
//     properly :/
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#ifndef __VECTOR_TEMPLATES_H__
#define __VECTOR_TEMPLATES_H__

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

template <typename T>
using matrix = std::vector<std::vector<T>>;

template <typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
	std::vector<T> result;
	for (size_t i = 0; i < a.size(); ++i) {
		result.push_back(a[i] + b[i]);
	}
	return result;
}

template <typename T>
std::vector<T> &operator+=(std::vector<T> &a, std::vector<T> const &b) {
	for (size_t i = 0; i < a.size(); ++i) {
		a[i] += b[i];
	}
	return a;
}

template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b) {
	std::vector<T> result;
	for (size_t i = 0; i < a.size(); ++i) {
		result.push_back(a[i] - b[i]);
	}
	return result;
}

template <typename T>
std::vector<T> &operator-=(std::vector<T> &a, std::vector<T> const &b) {
	for (size_t i = 0; i < a.size(); ++i) {
		a[i] -= b[i];
	}
	return a;
}

template <typename T>
std::vector<T> operator*(const std::vector<T> &a, const T b) {
	std::vector<T> result(a.begin(), a.end());
	for (auto &i : result) {
		i *= b;
	}
	return result;
}

template <typename T>
std::vector<T> operator*(const T b, const std::vector<T> &a) {
	return a * b;
}

template <typename T>
std::vector<T> &operator*=(std::vector<T> &a, T const &b) {
	for (auto &i : a) {
		i *= b;
	}
	return a;
}

template <typename T>
std::vector<T> operator%(std::vector<T> const &a, T const &b) {
	std::vector<T> result(a.begin(), a.end());
	for (auto &i : result) {
		i = (i + b) % b;
	}
	return result;
}

template <typename T>
std::vector<T> &operator%=(std::vector<T> &a, T const &b) {
	for (auto &i : a) {
		i = (i + b) % b;
	}
	return a;
}

template <typename T>
NTL::Vec<T> operator%(NTL::Vec<T> const &a, T const &b) {
	NTL::Vec<T> result;
	result.SetLength(a.length());
	for (long i = 0; i < a.length(); ++i) {
		result[i] = (a[i] + b) % b;
	}
	return result;
}

template <typename T>
NTL::Vec<T> &operator%=(NTL::Vec<T> &a, T const &b) {
	for (long i = 0; i < a.length(); ++i) {
		a[i] = (a[i] + b) % b;
	}
	return a;
}

// inner product
template <typename T>
T operator*(const std::vector<T> &a, const std::vector<T> &b) {
	return inner_product(a.cbegin(), a.cend(), b.cbegin(), (T)0);
}

#ifndef USE_FPLLL
template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
	std::string str;
	std::stringstream ss;
	ss << std::setprecision((int)os.precision()) << "[";
	for (auto &i : v) {
		ss << i << " ";
	}
	str = ss.str();
	str.erase(str.end() - 1, str.end());
	os << str << "]";
	return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os,
						 const std::vector<std::vector<T>> &v) {
	std::string str;
	std::stringstream ss;
	ss << std::setprecision((int)os.precision()) << "[";
	for (auto &i : v) {
		ss << i << std::endl;
	}
	str = ss.str();
	str.erase(str.end() - 1, str.end());
	os << str << "]";
	return os;
}

template <typename T>
std::istream &operator>>(std::istream &s, std::vector<T> &v) {
	T x;
	std::string token;
	char blank;

	s >> blank; // Gobble the first '['
	while (getline(s, token, ' ')) {
		std::istringstream input(token);
		input >> x;
		v.push_back(x);
	}
	// s >> blank; // Gobble the last ']'
	return s;
}

template <typename T>
std::istream &operator>>(std::istream &s, std::vector<std::vector<T>> &v) {
	std::vector<T> x;
	std::string token;
	char blank;

	s >> blank; // Gobble the first '['
	while (getline(s, token, '\n')) {
		std::istringstream input(token);
		input >> x;
		// s >> blank; //read ] after [...
		if (x.size() > 0) {
			v.push_back(x);
		}
		x.clear();
	}
	s >> blank; // Gobble the last ']'
	return s;
}
#endif // USE_FPLLL

#endif // __VECTOR_TEMPLATES_H__

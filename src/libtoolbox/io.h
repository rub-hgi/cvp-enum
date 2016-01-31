// ----------------------------------------------------------------------------
// Title      : IO Templates
// Project    : LWE Sampler
// ----------------------------------------------------------------------------
// File       : reduction.h
// Author     : Friedrich Wiemer <friedrich.wiemer@rub.de>
//              Elena Kirshanova <elena.kirshanova@rub.de>
// Company    : Ruhr-University Bochum
// Created    : 2015-10-28
// Last update: 2015-10-28
// ----------------------------------------------------------------------------
// Description:
//     Templated implementations for IO.
//     This needs to be in a header file, as c++ does not support templated
//     implementations properly :/
// ----------------------------------------------------------------------------
// Revisions  :
// Date        Version  Author  Description
// 2015-10-28  1.0      fwi     Created
// ----------------------------------------------------------------------------

#ifndef __IO_H__
#define __IO_H__

#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

/**
 * Write
 * \brief writes data to the file filename
 *
 * @params filename the file to write too
 * @params data to write
 */
template <typename T>
void Write(std::string const &filename, const T &data) {
	std::ofstream ofile;
	ofile.open(filename);
	if (ofile.is_open()) {
		ofile << std::setprecision(25) << data;
		ofile.close();
	} else {
		std::cout << "Unable to open file" << filename << std::endl;
	}
}

/**
 * Read
 * \brief reads data from the file filename and returns it
 *
 * @params filename the file to read the data from
 */
template <typename T>
T Read(std::string const &filename) {
	T data;
	std::ifstream ifile(filename);
	if (ifile.is_open()) {
		ifile >> data;
		ifile.close();
	} else {
		std::cout << "Unable to open file " << filename << std::endl;
	}
	return data;
}

#endif // __IO_H__

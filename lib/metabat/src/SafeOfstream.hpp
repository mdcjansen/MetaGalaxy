#ifndef UTILS
#define UTILS

#include <ios>
#include <iostream>
#include <string>
#include <stdio.h>

class SafeOfstream : public std::ofstream {
private:
	string realFileName, tmpFileName;
public:
	typedef std::ofstream Base;
	SafeOfstream( const char * filename, ios_base::openmode mode = ios_base::out): Base(filename, mode) {}
/*
	SafeOfstream( const char * filename, ios_base::openmode mode = ios_base::out) : realFileName(filename), tmpFileName(realFileName + ".tmp"), Base(tmpFileName.c_str(), mode) {
		std::cerr << "Opened " << tmpFileName << " with " << mode << std::endl;
	}
	virtual ~SafeOfstream() {
		this->~Base();
		rename(tmpFileName.c_str(), realFileName.c_str());
		std::cerr << "Renamed " << tmpFileName << " to " << realFileName << std::endl;
	}
*/
};


#endif

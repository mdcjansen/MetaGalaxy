/*
 * OpenMP.h
 *
 *  Generic Header to include basic OpenMP functions
 *  when OpenMP is not enabled
 *
 */

#ifndef OPENMP_H_
#define OPENMP_H_


#ifdef _OPENMP

#include <omp.h>

#else

#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_in_parallel() 0

void omp_set_num_threads(int i) {}
void omp_set_nested(int i) {}

#endif

#endif /* OPENMP_H_ */

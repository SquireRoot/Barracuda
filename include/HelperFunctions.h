#ifndef HELPER_FUNCIONS_H
#define HELPER_FUNCIONS_H

#include <stdio.h>
#include "Config.h"

// index for 3D
#define I3D(x,y,z,nx,ny) 	((nx)*((z)*(ny) + (y)) + x)
#define I3Df(i3d,i,ni) 		(((i3d)*(ni)) + i)

// printing helpers to view arrays
void printRealArray(real *arr, int sizeX, int sizeY, int sizeZ);
void printIntArray(int *arr, int sizeX, int sizeY, int sizeZ);
void printVectArray(vect *arr, int sizeX, int sizeY, int sizeZ);
void printFArray(real *arr, int sizeI, int sizeX, int sizeY, int sizeZ);

// Cuda error checking functions
#define cudaSafeCall(err) __cudaSafeCall(err, __FILE__, __LINE__)
#define cudaCheckError() __cudaCheckError(__FILE__, __LINE__)
 
inline void __cudaSafeCall(cudaError err, const char *file, const int line) {
	if ( cudaSuccess != err ) {
		fprintf(stderr, "cudaSafeCall() failed at %s:%i : %s\n",
				file, line, cudaGetErrorString(err));
		exit(-1);
	}
	return;
}
 
inline void __cudaCheckError(const char *file, const int line) {
	cudaError err = cudaGetLastError();
	if (cudaSuccess != err) {
		fprintf(stderr, "cudaCheckError() failed at %s:%i : %s\n",
				file, line, cudaGetErrorString(err));
		exit(-1);
	}
	 
	// More careful checking. However, this will affect performance.
	// Comment away if needed.
	err = cudaDeviceSynchronize();
	if(cudaSuccess != err) {
		fprintf(stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
				file, line, cudaGetErrorString(err));
		exit(-1);
	}
	return;
}

#endif
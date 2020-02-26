#ifndef KERNELS_CUH
#define KERNELS_CUH

#include "Config.h"
#include "HelperFunctions.h"

__global__ void lbmKernel(LatticeParams* params, ObservableArrays* arrays,
						   real* fNew, real* fOld, int* nodeType);
__global__ void initKernel(LatticeParams* params, ObservableArrays* arrays,
						   real* fNew, real* fOld, int* nodeType);
__device__ real dot(vect a, vect b);
__device__ real fEQ(real w, real rho, vect e, vect u);

#endif
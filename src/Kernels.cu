#include "Kernels.cuh"
#include <stdio.h>

__global__ void lbmKernel(LatticeParams* params, ObservableArrays* arrays,
						  real* fNew, real*fOld, int* nodeType) {
	int x = threadIdx.x;
	int y = blockIdx.x;
	int z = blockIdx.y;
	int index = I3D(x, y, z, SIZE_X, SIZE_Y);

	real f[NUM_DIRECTIONS];

	/* --- Streaming: pull from 0 to 1 --- */
	// i = 0 case is easy
	f[0] = fOld[index*NUM_DIRECTIONS];
	// loop through remaining distributions
	for (int i = 1; i < NUM_DIRECTIONS; i++) {
		// default periodic boundary in each direction
		int indexei = I3D(((x - int(params->ei[i].x) + SIZE_X)%SIZE_X),
						  ((y - int(params->ei[i].y) + SIZE_Y)%SIZE_Y),
						  ((z - int(params->ei[i].z) + SIZE_Z)%SIZE_Z),
						  SIZE_X, SIZE_Y);

		if (nodeType[index] == 0) {
			switch (nodeType[indexei]) {
			case 0:
				f[i] = fOld[I3Df(indexei, i, NUM_DIRECTIONS)];
				break;
			case 1:
				f[i] = fOld[I3Df(index, params->bbi[i], NUM_DIRECTIONS)];
				break;
			default:
				f[i] = fOld[I3Df(index, i, NUM_DIRECTIONS)];
			}

		} else {
			f[i] = fOld[I3Df(index, i, NUM_DIRECTIONS)];
		}
	}

	real rho = 0.0;
	for (int i = 0; i < NUM_DIRECTIONS; i++) {
		rho = rho + f[i];
	}
	arrays->rho[index] = rho;

	vect u = {0.0, 0.0, 0.0};
	for (int i = 0; i < NUM_DIRECTIONS; i++) {
		u.x = u.x + (f[i]*params->ei[i].x)/rho;
		u.y = u.y + (f[i]*params->ei[i].y)/rho;
		u.z = u.z + (f[i]*params->ei[i].z)/rho;
	}
	arrays->u[index] = u;

	for (int i = 0; i < NUM_DIRECTIONS; i++) {
		real equilibrium = fEQ(params->wi[i], rho,
							   params->ei[i], u);
		fNew[I3Df(index, i, NUM_DIRECTIONS)] = f[i] - (f[i] - equilibrium)/TAU;
		// fNew[I3Df(index, i, NUM_DIRECTIONS)] = f[i];
	}
}


__global__ void initKernel(LatticeParams* params, ObservableArrays* arrays, 
						   real* fNew, real* fOld, int* nodeType) {
	int x = threadIdx.x;
	int y = blockIdx.x;
	int z = blockIdx.y;

	int index = I3D(x, y, z, SIZE_X, SIZE_Y);

	// intialize boundary
	if (y == 0 || y == SIZE_Y - 1) {
		nodeType[index] = 1;
	} else {
		nodeType[index] = 0;
	}

	// initialize u
	arrays->u[index].x = INIT_U_X;
	arrays->u[index].y = INIT_U_Y;
	arrays->u[index].z = INIT_U_Z;

	// initialize rho
	arrays->rho[index] = INIT_RHO;

	// initialize distributions
	for (int i = 0; i < NUM_DIRECTIONS; i++) {
		real eqDist = fEQ(params->wi[i], INIT_RHO,
						  params->ei[i], arrays->u[index]);
		fOld[I3Df(index, i, NUM_DIRECTIONS)] = eqDist;
		fNew[I3Df(index, i, NUM_DIRECTIONS)] = eqDist;
	}
}

__device__ real fEQ(real w, real rho, vect e, vect u){
	real eDotu = dot(e,u);
	real eqDist = rho*w*(1.0 + 3.0*eDotu + 4.5*eDotu*eDotu - 1.5*dot(u,u));
	return eqDist;
}

__device__ real dot(vect a, vect b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
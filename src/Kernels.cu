#include "Kernels.cuh"
#include <stdio.h>

__global__ void lbmKernel(LatticeParams* params, ObservableArrays* arrays,
						  real* fNew, real*fOld, int* nodeType) {
	int x = threadIdx.x;
	int y = blockIdx.x;
	int z = blockIdx.y;
	int index = I3D(x, y, z, SIZE_X, SIZE_Y);

	if (nodeType[index] == 0) {

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

		}


		real rho = 0.0;
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			rho = rho + f[i];
		}

		vect u = {0.0, 0.0, 0.0};
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			u.x = u.x + (f[i]*params->ei[i].x)/rho;
			u.y = u.y + (f[i]*params->ei[i].y)/rho;
			u.z = u.z + (f[i]*params->ei[i].z)/rho;
		}

		#ifdef ENABLE_FORCE
		const vect force = {FORCE_X, FORCE_Y, FORCE_Z};
		real* fp = (real*)(&force);
		real* up = (real*)(&u);

		real source[NUM_DIRECTIONS];
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			real* cip = (real*)(&(params->ei[i]));
			source[i] = 0.0;
			for (int a = 0; a < 3; a++) {
				for (int b = 0; b < 3; b++) {
					source[i] = source[i] + 
					fp[a]*(3.0*cip[a] + 9.0*up[b]*(cip[a]*cip[b] - (a == b ? 1.0 : 0.0)/3.0));
				}
			}
			source[i] = source[i]*params->wi[i];

			rho = rho + source[i]/2;

			real calpha = 0.0;
			for (int a = 0; a < 3; a++) {
				calpha = calpha + cip[a];
			}

			u.x = u.x + source[i]*calpha*0.5/rho;
			u.y = u.y + source[i]*calpha*0.5/rho;
			u.z = u.z + source[i]*calpha*0.5/rho;
		}

		#define TAU_BAR  (TAU + 0.5)
		#else
		#define TAU_BAR  (TAU)
		#endif

		arrays->rho[index] = rho;
		arrays->u[index] = u;

		#ifdef SINGLE_RELAXATION
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			real equilibrium = fEQ(params->wi[i], rho,
								   params->ei[i], u);
			f[i] = f[i] - (f[i] - equilibrium)/TAU_BAR;
		}
		#endif

		#ifdef ENABLE_FORCE
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			f[i] = f[i] + (1 - 1/(2*TAU_BAR))*source[i];
		}
		#endif

		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			fNew[I3Df(index, i, NUM_DIRECTIONS)] = f[i];
		}
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
	}

	#ifdef ENABLE_CYLINDER
	else if ((x - CYLINDER_X0)*(x - CYLINDER_X0) +
			 (y - CYLINDER_Y0)*(y - CYLINDER_Y0) 
			  < CYLINDER_R*CYLINDER_R) { 
		nodeType[index] = 1;
	}
	#endif

	else {
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
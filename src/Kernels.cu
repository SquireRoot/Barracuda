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
			
			#ifdef PERIODIC_BOUNDARY
			int indexei = I3D(((x - int(params->ei[i].x) + SIZE_X)%SIZE_X),
							  ((y - int(params->ei[i].y) + SIZE_Y)%SIZE_Y),
							  ((z - int(params->ei[i].z) + SIZE_Z)%SIZE_Z),
							  SIZE_X, SIZE_Y);
			#else
			int indexei = I3D((x - int(params->ei[i].x)),
							  (y - int(params->ei[i].y)),
							  (z - int(params->ei[i].z)),
							  SIZE_X, SIZE_Y);
			#endif

			switch (nodeType[indexei]) {
			case 0: // fluid
				f[i] = fOld[I3Df(indexei, i, NUM_DIRECTIONS)];
				break;
			case 1: // bounce back
				f[i] = fOld[I3Df(index, params->bbi[i], NUM_DIRECTIONS)];
				break;
			case 2: // slipping bounce back
				f[i] = fOld[I3Df(index, params->sbbxi[i], NUM_DIRECTIONS)];				
				break;
			case 3: // inlet
			{ // declare a local scope
				// This code assumes there are no wet inlet nodes on the edge of the domain
				const vect xDir = {1.0, 0.0, 0.0};
				const vect yDir = {0.0, 1.0, 0.0};
				// vect zDir = {0.0, 0.0, 1.0};

				// 1.0 if first order accurate asymetric, 0.5 if second order accurate symetric
				real xDerivativeFactor = (fabsf(dot(xDir, params->ei[i])) > 0 ? 1.0 : 0.5);
				real yDerivativeFactor = (fabsf(dot(yDir, params->ei[i])) > 0 ? 1.0 : 0.5);
				// real zDerivativeFactor = (fabsf(dot(zDir, params->ei[i])) > 0 ? 1.0 : 0.5);

				real gradRhoU = (-arrays->u[I3D(x - 1, y, z, SIZE_X, SIZE_Y)].x 
								 +arrays->u[I3D(x + 1, y, z, SIZE_X, SIZE_Y)].x)*xDerivativeFactor;
				gradRhoU	 += (-arrays->u[I3D(x, y - 1, z, SIZE_X, SIZE_Y)].y 
								 +arrays->u[I3D(x, y + 1, z, SIZE_X, SIZE_Y)].y)*yDerivativeFactor;
				// gradRhoU	 += (-arrays->u[I3D(x, y, z - 1, SIZE_X, SIZE_Y)].z 
				// 				 +arrays->u[I3D(x, y, z + 1, SIZE_X, SIZE_Y)].z)*zDerivativeFactor;
				gradRhoU *= arrays->rho[index];

				// Set i = 0 cause its easy
				fNew[index*NUM_DIRECTIONS] =
						fEQ(params->wi[0], INIT_RHO, params->ei[0], arrays->u[index])
						+TAU*params->wi[0]*gradRhoU;
				// Compute remaining distributions
				for (int gi = 1; gi < NUM_DIRECTIONS; gi++) {
					real eiDotGradEiDotRhoU = params->ei[gi].x
							*(-dot(arrays->u[I3D(x - 1, y, z, SIZE_X, SIZE_Y)], params->ei[gi]) 
							  +dot(arrays->u[I3D(x + 1, y, z, SIZE_X, SIZE_Y)], params->ei[gi]))
							*xDerivativeFactor;
					eiDotGradEiDotRhoU += params->ei[gi].y
							*(-dot(arrays->u[I3D(x, y - 1, z, SIZE_X, SIZE_Y)], params->ei[gi]) 
							  +dot(arrays->u[I3D(x, y + 1, z, SIZE_X, SIZE_Y)], params->ei[gi]))
							*yDerivativeFactor;
					// eiDotGradEiDotRhoU += params->ei[gi].z
					// 		*(-dot(arrays->u[I3D(x, y, z - 1, SIZE_X, SIZE_Y)], params->ei[gi]) 
					// 		  +dot(arrays->u[I3D(x, y, z + 1, SIZE_X, SIZE_Y)], params->ei[gi]))
					// 		*zDerivativeFactor;
					eiDotGradEiDotRhoU *= arrays->rho[index];


					fNew[I3Df(index, gi, NUM_DIRECTIONS)] =
							fEQ(params->wi[gi], INIT_RHO, params->ei[gi], arrays->u[index])
							-TAU*params->wi[gi]*(3.0*eiDotGradEiDotRhoU - gradRhoU);
				}

			}
				return;

			case 4: {// outlet

				int indexl = I3D(x-1, y, z, SIZE_X, SIZE_Y);
				f[i] = fOld[I3Df(indexl, i, NUM_DIRECTIONS)];
			}	break;

			default:
				f[i] = fOld[I3Df(index, i, NUM_DIRECTIONS)];
			}

		} // end local scope declaration


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
		#endif

		arrays->rho[index] = rho;
		arrays->u[index] = u;

		#ifdef SINGLE_RELAXATION
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			real equilibrium = fEQ(params->wi[i], rho, params->ei[i], u);
			f[i] = f[i] - (f[i] - equilibrium)/TAU;
		}
		#elif defined TWO_RELAXATION
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			real fp = 0.5*(f[i] + f[params->bbi[i]]);
			real fm = 0.5*(f[i] - f[params->bbi[i]]);
			real feqp = 0.5*(fEQ(params->wi[i], rho, params->ei[i], u)
						     +fEQ(params->wi[params->bbi[i]],
						          rho, params->ei[params->bbi[i]], u));
			real feqm = 0.5*(fEQ(params->wi[i], rho, params->ei[i], u)
						     -fEQ(params->wi[params->bbi[i]],
						          rho, params->ei[params->bbi[i]], u));
			f[i] = ((1.0 - OMEGAP)*fp + OMEGAP*feqp)
				  +((1.0 - OMEGAM)*fm + OMEGAM*feqm);
			//f[i] = fp + fm - OMEGAP*(fp + fm - (feqp + feqm));
		}
		#endif

		#ifdef ENABLE_FORCE
		for (int i = 0; i < NUM_DIRECTIONS; i++) {
			f[i] = f[i] + (1 - 1/(2*TAU))*source[i];
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
		nodeType[index] = 2;
	}

	#ifdef ENABLE_CYLINDER
	else if ((x - CYLINDER_X0*SIZE_X)*(x - CYLINDER_X0*SIZE_X) +
			 (y - CYLINDER_Y0*SIZE_Y)*(y - CYLINDER_Y0*SIZE_Y) 
			  < CYLINDER_R*CYLINDER_R) { 
		nodeType[index] = 1;
	}
	#endif

	else {
		nodeType[index] = 0;
	}

	#ifdef ENABLE_INLET
	if (x == 0) {
		nodeType[index] = 3;
	} 
	if (x == SIZE_X - 1) {
		nodeType[index] = 3;
	}
	#endif

	/* Initialize u */
	if (nodeType[index] != 1) {
		arrays->u[index] = {INIT_U_X, INIT_U_Y, INIT_U_Z};
	} else {
		arrays->u[index] = {0.0, 0.0, 0.0};
	}

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
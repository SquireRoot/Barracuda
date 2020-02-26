#include <stdio.h>
#include "Config.h"
#include "HelperFunctions.h"
#include "Kernels.cuh"

LatticeParams* hParams;
LatticeParams* dParams;

ObservableArrays* hArrays;
ObservableArrays* dArrays;

real* dFNew;
real* dFOld;
int* nodeType;

real* hFCurrent;

int main(int argc, char **argv) {
	printf("Barracuda v%2.2f\n", VERSION);
	printf("Georgia Institute of Technology Department of Physics - CHAOS Lab\n");
	printf("Author: Evan Newman - evanmnewman1@gmail.com\n");
	printf("Special Help from Abouzar Kaboudian and Flavio Fenton\n");

	/* Write the config file */
	FILE *fp = fopen("./DATA/config.txt", "w");
	fprintf(fp, "SizeX = %d\n", SIZE_X);
	fprintf(fp, "SizeY = %d\n", SIZE_Y);
	fprintf(fp, "SizeZ = %d\n", SIZE_Z);
	fprintf(fp, "Iterations = %d\n", NUM_TIMESTEPS);
	fprintf(fp, "OutputInterval = %d\n", OUTPUT_INTERVAL);
	#if OUTPUT_START < OUTPUT_INTERVAL
		fprintf(fp, "OutputBegin = %d\n", OUTPUT_INTERVAL);
	#else
		fprintf(fp, "OutputBegin = %d\n", OUTPUT_START);
	#endif
	fclose(fp);

	/* Allocate memory on host and device*/
	hParams = new LatticeParams;
	hArrays = new ObservableArrays;

	cudaSafeCall(cudaMalloc((void**)&dParams, sizeof(LatticeParams)));
	cudaSafeCall(cudaMemcpy(dParams, hParams, sizeof(LatticeParams),
							cudaMemcpyHostToDevice));

	cudaSafeCall(cudaMalloc((void**)&dArrays, sizeof(ObservableArrays)));
	cudaSafeCall(cudaMemcpy(dArrays, hArrays, sizeof(ObservableArrays),
							cudaMemcpyHostToDevice));

	cudaSafeCall(cudaMalloc((void**)&dFNew, sizeof(real)*DOMAIN_SIZE*NUM_DIRECTIONS));
	cudaSafeCall(cudaMalloc((void**)&dFOld, sizeof(real)*DOMAIN_SIZE*NUM_DIRECTIONS));
	cudaSafeCall(cudaMalloc((void**)&nodeType, sizeof(int)*DOMAIN_SIZE));	
	hFCurrent = (real*)malloc(sizeof(real)*DOMAIN_SIZE*NUM_DIRECTIONS);

	/* Initialize arrays on host */
	dim3 grid = dim3(SIZE_Y, SIZE_Z, 1);
	dim3 block = dim3(SIZE_X, 1, 1);
	initKernel<<<grid, block>>>(dParams, dArrays, dFNew, dFOld, nodeType);
	cudaDeviceSynchronize();
	cudaCheckError();

	printf("Finished initialization\n");
	
	/* Loop */
	for (int i = 1; i <= NUM_TIMESTEPS; i++) {

		lbmKernel<<<grid, block>>>(dParams, dArrays, dFNew, dFOld, nodeType);
		cudaDeviceSynchronize();
		cudaCheckError();
		
		if (i >= OUTPUT_START && i%OUTPUT_INTERVAL == 0) { // output data 
			cudaSafeCall(cudaMemcpy(hArrays, dArrays, sizeof(ObservableArrays),
									cudaMemcpyDeviceToHost));

			FILE *fpv;
		    char filename[sizeof "./DATA/curl/file3000000.dat"];
		    sprintf(filename, "./DATA/v/file%07d.dat", i);
		    fpv = fopen(filename, "w+b");
		    fwrite(hArrays->u, sizeof(vect), DOMAIN_SIZE, fpv);
		    fclose(fpv);

			printf("Wrote data at t = %i\n", i);
		}

		real* temp = dFOld;
		dFOld = dFNew;
		dFNew = temp;
	}

	return 0;
}
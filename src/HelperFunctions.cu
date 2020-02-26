#include "HelperFunctions.h"

void printRealArray(real *arr, int sizeX, int sizeY, int sizeZ) {
	for (int z = 0; z < sizeZ; z++) {
		printf("z = %i:\n", z);
		for (int y = sizeY - 1; y >= 0; y--) {
			printf("y = %2i |  ", y);
			for (int x = 0; x < sizeX; x++) {
				printf("%2.3f ", arr[I3D(x, y, z, sizeX, sizeY)]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void printIntArray(int *arr, int sizeX, int sizeY, int sizeZ) {
	for (int z = 0; z < sizeZ; z++) {
		printf("z = %i:\n", z);
		for (int y = sizeY - 1; y >= 0; y--) {
			printf("y = %2i |  ", y);
			for (int x = 0; x < sizeX; x++) {
				printf("%4i ", arr[I3D(x, y, z, sizeX, sizeY)]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void printVectArray(vect *arr, int sizeX, int sizeY, int sizeZ) {
	for (int z = 0; z < sizeZ; z++) {
		printf("z = %i:\n", z);
		for (int y = sizeY - 1; y >= 0; y--) {
			printf("y = %2i:\n", y);
			
			printf("        x |  ");
			for (int x = 0; x < sizeX; x++) {
				printf("%2.3f ", arr[I3D(x, y, z, sizeX, sizeY)].x);
			}
			printf("\n");

			printf("        y |  ");
			for (int x = 0; x < sizeX; x++) {
				printf("%2.3f ", arr[I3D(x, y, z, sizeX, sizeY)].y);
			}
			printf("\n");

			printf("        z |  ");
			for (int x = 0; x < sizeX; x++) {
				printf("%2.3f ", arr[I3D(x, y, z, sizeX, sizeY)].z);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void printFArray(real *arr, int sizeI, int sizeX, int sizeY, int sizeZ) {
	for (int z = 0; z < sizeZ; z++) {
		printf("z = %i:\n", z);
		for (int y = sizeY - 1; y >= 0; y--) {
			printf("y = %2i:\n", y);
			for (int i = 0; i < sizeI; i++) {
				printf("        i = %2i |  ", i);
				for (int x = 0; x < sizeX; x++) {
					printf("%2.3f ", arr[I3Df(I3D(x, y, z, sizeX, sizeY), i, sizeI)]);
				}
				printf("\n");
			}
		}
		printf("\n");
	}
}
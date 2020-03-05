#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include "UserConfig.h"

#define VERSION (1.0)

#ifdef DOUBLE_PRECISION
	typedef double real;
	typedef double3 vect;
#else
	typedef float real;
	typedef float3 vect;
#endif

#ifdef ENABLE_FORCE
	#define TAU (3.0*VISCOSITY + 1.0)
#else
	#define TAU ((3.0*VISCOSITY) + 0.5)
#endif

#ifdef TWO_RELAXATION
	#define OMEGAP (1.0/TAU)
	#define OMEGAM (1.0/((MAGIC_PARAMETER)/(TAU - 0.5) + 0.5)) 
#endif

#ifdef D2_Q9
	/* index directions
	6   2   5
	  \ | /
	3 - 0 - 1
	  / | \
	7   4   8
	*/
	#define NUM_DIRECTIONS (9)
	// CUDA C allows us to have default initial values in a struct and use new
	typedef struct LatticeParamsStruct {
		vect ei[9] = {{ 0.0,  0.0,  0.0}, // 0
					  { 1.0,  0.0,  0.0}, // 1
					  { 0.0,  1.0,  0.0}, // 2
					  {-1.0,  0.0,  0.0}, // 3
					  { 0.0, -1.0,  0.0}, // 4
					  { 1.0,  1.0,  0.0}, // 5
					  {-1.0,  1.0,  0.0}, // 6
					  {-1.0, -1.0,  0.0}, // 7
					  { 1.0, -1.0,  0.0}  // 8
		};

		int bbi[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
		int sbbxi[9] = {0, 3, 4, 1, 2, 8, 7, 6, 5};

		real wi[9] = {4.0/9.0,
					  1.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,
					  1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
	} LatticeParams;
#elif defined D3_Q15	
	// d3_q15 definitions...
#endif

#define DOMAIN_SIZE (SIZE_X*SIZE_Y*SIZE_Z)

typedef struct ObservableArraysStruct {
	vect u[DOMAIN_SIZE];
//	vect dp[DOMAIN_SIZE];
	real rho[DOMAIN_SIZE];
//	vect curl[DOMAIN_SIZE];
} ObservableArrays;

#endif

// Lex and Yacc - tokenizer and grammer
// stress free boundary conditions 
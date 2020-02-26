#ifndef USER_CONFIG_H
#define USER_CONFIG_H

/* Control Parameters */
#define NUM_TIMESTEPS    (1000)
#define OUTPUT_INTERVAL  (10)
#define OUTPUT_START     (1)

/* Doman Parameters */
#define SIZE_X    (512)
#define SIZE_Y    (512)
#define SIZE_Z    (1)

/* Lattice type */
#define D2_Q9
//#define D3_Q15

/* Precision - default is single precision*/
//#define DOUBLE_PRECISION

// STL Units per lattice spacing
#define STL_TO_LATTICE (1.0)

// Initial Parameters
#define INIT_U_X   (0.01)
#define INIT_U_Y   (0.0)
#define INIT_U_Z   (0.0)

#define INIT_RHO   (1.0)

#define VISCOSITY  (0.001)

#endif
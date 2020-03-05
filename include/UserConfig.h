#ifndef USER_CONFIG_H
#define USER_CONFIG_H

/* Control Parameters */
#define NUM_TIMESTEPS    (500000)
#define OUTPUT_INTERVAL  (500)
#define OUTPUT_START     (1)

/* Doman Parameters */
#define SIZE_X    (512)
#define SIZE_Y    (512)
#define SIZE_Z    (1)

/* Lattice type */
#define D2_Q9
//#define D3_Q15

/* Relaxation Scheme */
#define SINGLE_RELAXATION

//#define TWO_RELAXATION
//#define MAGIC_PARAMETER   (0.25)

/* Precision - default is single precision*/
//#define DOUBLE_PRECISION

/* Initial Velocity */
#define INIT_U_X   (0.02)
#define INIT_U_Y   (0.0)
#define INIT_U_Z   (0.0)

/* Initial Density */
#define INIT_RHO   (1.0)

/* Kinematic Viscosity */
#define VISCOSITY  (0.0003)

/* Constant Force */
//#define ENABLE_FORCE
#define FORCE_X    (0.00000002)
#define FORCE_Y    (0.0)
#define FORCE_Z    (0.0)

/* Periodic Boundary Conditions */
//#define PERIODIC_BOUNDARY

/* Inlet Boundary at x = 0, velocity willbe INIT_U */
#define ENABLE_INLET

/* Cylinder Boundary */
#define ENABLE_CYLINDER
#define CYLINDER_X0   (0.25)
#define CYLINDER_Y0   (0.5)
#define CYLINDER_R    (13)


#endif


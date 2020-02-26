#ifndef USER_CONFIG_H
#define USER_CONFIG_H

/* Control Parameters */
#define NUM_TIMESTEPS    (500000)
#define OUTPUT_INTERVAL  (500)
#define OUTPUT_START     (50000)

/* Doman Parameters */
#define SIZE_X    (512)
#define SIZE_Y    (512)
#define SIZE_Z    (1)

/* Lattice type */
#define D2_Q9
//#define D3_Q15

/* Relaxation Scheme */
#define SINGLE_RELAXATION

/* Precision - default is single precision*/
//#define DOUBLE_PRECISION

/* Initial Velocity */
#define INIT_U_X   (0.0)
#define INIT_U_Y   (0.0)
#define INIT_U_Z   (0.0)

/* Initial Density */
#define INIT_RHO   (1.0)

/* Kinematic Viscosity */
#define VISCOSITY  (0.001)

/* Constant Force */
#define ENABLE_FORCE
#define FORCE_X    (0.00000002)
#define FORCE_Y    (0.0)
#define FORCE_Z    (0.0)

/* Cylinder Boundary */
#define ENABLE_CYLINDER
#define CYLINDER_X0   (128)
#define CYLINDER_Y0   (256)
#define CYLINDER_R   (5)


#endif
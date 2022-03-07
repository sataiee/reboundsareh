/* This file contains the  global variables for reboundsareh */
#define MAX_SIZE 50
#define THICKNESSSMOOTHING 0.6
#define ADIABATICINDEX     1.0 //for a locally isothermal it is 1.0, must be changed for a non-isothermal model
#define DIFFUSIVITY        0.0 //This will be used for discs with thermal diffusion
#define MAX1D              2000

double Gamma[MAX_SIZE], GammaL[MAX_SIZE], GammaHSB[MAX_SIZE], GammaHSE[MAX_SIZE], Gamma0[MAX_SIZE];
double GammaLB[MAX_SIZE], GammaLE[MAX_SIZE], GammaCTotal[MAX_SIZE], sg_wrong_factor[MAX_SIZE]; 
double GammaStar[MAX_SIZE];
double tau_ecc[MAX_SIZE], tau_mig[MAX_SIZE], tau_a[MAX_SIZE];
init_planet *planets;
double tstart;

double alpha_sigma, sigma0, h0, flaring, tmax, tdump, alpha_visc;
int without_sg_term = 0, without_e_damping = 0, isothermal = 0, unsaturated = 0, torque_to_acc = 0;
/*
   * without_sg_term = 1 if we want to exclude the sg term (as in hydro) otherwise it should be always 0
   * without_e_damping = 1 if we want to exclude the eccentricity damping otherwise it should be always 0
   * isothermal = 0 for locally isothermal models, and non-zero for radiative ones
   * torque_to_acc = 0 means the method of converting torque to acceleration is Hano (default rebound)
   * torque_to_acc = 1 means the method of converting torque to acceleration is Willy (as Willy's calculation)
   * torque_to_acc = 2 means the method of converting torque to acceleration is John (as Papaloizou 2011)
*/

/* These store the radii and surface density array when the surface density is read 
 * from a file named sigma.dat. It must be an ascii file with two column. The first column
 * should be r and the second sigma
*/
double sigma_array[MAX1D], rad_array[MAX1D];
int read_sig_from_file = 1; // 0: read from file, 1:default power-low

int stellar_torque;  //if 0, stellar torque is included
double R_star_to_sun;  //Star-to-sun radius ratio
double R_co;           //corotation radius in au

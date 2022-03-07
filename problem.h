/* --------type definitions------- */
typedef struct basic_variable {
  char *name;
  double value;
  int read;  // 0 if the variable hasn't read and 1 if it has
} basic_variable;

typedef struct init_planet {
  double mass;			// masses 
  double a;			        // initial semi-major axis
  double e;        // eccentricity
  char    mig_type[20];  //migration type
  double tau_a;  //tau_a
  double tau_e;  //tau_e
  double inc;  //inclination
  double Omega; //longitude of the ascending node
  double omega; //argument of the pericentre
  double f;
  double t_pert; //if not 0, a perturbation is applied on the planet's velocities 
  double a_pert; // amplitude of the perturbation
} init_planet;

/* --------function declerations ---------*/
void additional_forces(struct reb_simulation* const r);
void heartbeat(struct reb_simulation* const r);
void read_variables();
void accel_from_formula(struct reb_simulation* const r, double alpha_sigma, double beta_temperature);
double calculate_aspect_ratio(double dist);
double calculate_sigma_at_rp(double dist); 
double calculate_sigma_slope(double dist); 
double calculate_temperature_slope(double dist);
double caclulate_viscosity(double dist, double h);
double Ffunc(double p);
double Gfunc(double p);
double Kfunc(double p);
void read_planets(init_planet *planets, char *which_input, int nplanets);
int find_planets_number(char *which_input);
void read_sigma();


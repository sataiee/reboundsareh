/* @author: Sareh Ataiee   *
 * This setup, can perfrom 
 * 1) Migration of a planet in a disc with given properties
 *    using Paardekooper+fendyke models and given surface density profile
 * 2) Migration of the planets at the disc inner edge
 *    using different profiles
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <time.h>
#include <dirent.h>
#include "rebound.h"
#include "problem.h"
#include "global.h"

int counter_archive = 0;
double tmax  = 2e6*2*M_PI;
double tdump = 10*2*M_PI;

int main(int argc, char* argv[]){
    srand((long)time(NULL)); /* initialize rand() */
    int nplanets, i;
    char read_from[20] = "restartdisc"; //can be fargo, rebound, disc, or restartdisc
    if (strcmp(read_from, "fargo") == 0){
        printf("It is pure n-body.\n");
        printf("Reading the planets orbital elements from fargo outputs...\n");
    } else if (strcmp(read_from, "rebound") == 0){
        printf("It is pure n-body.\n");
        printf("Reading the planets orbital elements from rebound outputs...\n");
    } else if (strcmp(read_from, "disc") == 0){
        printf("This run includes migration.\n");
        printf("Reading the planets properties from planets.par file...\n");
    } else if (strcmp(read_from, "restartdisc") == 0){
        printf("This run is a restart from a previous migration model.\n");
    } else {
        printf("Invalid value for read_from :-/\n");
        exit(1);
    }
    
    nplanets = find_planets_number(read_from);  //number of planets
    planets = (init_planet *)malloc(sizeof(init_planet) * nplanets);
    read_planets(planets, read_from, nplanets);

    // Initializing the rebound structure
    struct reb_simulation* r = reb_create_simulation();
    
    if ((strcmp(read_from,"disc") == 0) || (strcmp(read_from,"fargo") == 0)){
        // Delete previous output
        system("rm -v *.txt");
        system("rm -v *.bin");
    }
    
    // Setup constants
    r->dt             = 1e-3*2*M_PI;        // initial timestep.
    //r->integrator     = REB_INTEGRATOR_IAS15;
    r->integrator     = REB_INTEGRATOR_WHFAST;

    if ((strcmp(read_from, "restartdisc") == 0) || strcmp(read_from, "disc") == 0){
        int n_vars=14;
        // reading and setting the variables
        read_variables(n_vars);

        //If read_sig_from_file, here we fill the rad and sigma array
        if (read_sig_from_file == 0)
            read_sigma();

        // Setup callback function for velocity dependent forces.
        r->additional_forces     = additional_forces;
        r->force_is_velocity_dependent = 1;
    }

    // Setup callback function for outputs.
    r->heartbeat        = heartbeat;

    struct reb_particle star = {0}; //star
    star.m = 1.0; 
    // Adding the star to the simulation
    reb_add(r, star);
    
    for(i=0; i<nplanets; i++) {
        struct reb_particle p = {0}; // setting all parameters of the planet to zero
    	p = reb_tools_orbit_to_particle(r->G, star, planets[i].mass, planets[i].a, planets[i].e, planets[i].inc, \
        planets[i].Omega, planets[i].omega, planets[i].f);
        reb_add(r, p);
    }
    
    //Move the calculation to the center of mass
    reb_move_to_com(r);
    
    // Do the integration
    reb_integrate(r, tmax);
    
    printf("\n");
    
    //free the planets array
    free(planets);

    //Make a binary simulation archive 
    reb_output_binary(r, "restart_final.bin");

    // Write the final orbital elements
    //reb_output_orbits(r, "final_element.txt");

    //return the centre of mass 
    reb_get_com(r);

    exit(0);
}



void read_variables(int n){
    char line[350], variable[100];
    double value;
    char *names[14] = {"alpha_sigma", "sigma0", "h0", "flaring", "without_sg_term", \
                        "without_e_damping", "alpha_visc", "isothermal", "unsaturated",\
                        "torque_to_acc", "read_sig_from_file", "stellar_torque", "R_star_to_sun",\
                        "R_co"};
    int i, check_n_vars=0;   
    
    // array for reading the disc and planet parameters
    basic_variable vars[n];
    
    // Initializing variable array
    for (i=0; i<n; i++){
        vars[i].name = names[i];
        vars[i].value = 0.0;
        vars[i].read = 0;  
    }
        
    //Read the parameter file
    FILE *input = fopen("parameters.par", "r");
    if (input == NULL)
        reb_exit("The parameter file does not exist.\n");
    
    // Read the variables
    while(fgets(line, 349, input) != NULL){
        if ((line[0] != '#') && (line[0] != '\n')){
            check_n_vars++;
            sscanf(line, "%s %lf", variable, &value);
            for (i=0; i<n; i++){
                if (strcmp(vars[i].name,variable) == 0){
                    vars[i].value = value;
                    if (vars[i].read == 0.0)
                        vars[i].read = 1.0;
                    else{
                        printf("%s, %e, %d\n", vars[i].name, vars[i].value, vars[i].read);
                        reb_exit("One variable is repeated more than once.\n");
                    }
                }
            }
        }
    }
    if (check_n_vars != n){
        printf("read variables=%d, number of variables must be %d\n", check_n_vars, n);
        reb_exit("Some variables are missing or extra in the parameters file\n");
    }
    fclose(input);
    
    for (int i=0; i<n; i++){
        if (strcmp(vars[i].name, "alpha_sigma")==0)     alpha_sigma = vars[i].value;
        if (strcmp(vars[i].name, "sigma0")==0)          sigma0 = vars[i].value;
        if (strcmp(vars[i].name, "h0")==0)              h0 = vars[i].value;
        if (strcmp(vars[i].name, "flaring")==0)         flaring = vars[i].value;
        if (strcmp(vars[i].name, "without_sg_term")==0) without_sg_term = (int) vars[i].value;
        if (strcmp(vars[i].name, "without_e_damping")==0) without_e_damping = (int) vars[i].value;
        if (strcmp(vars[i].name, "alpha_visc")==0)      alpha_visc = vars[i].value;
        if (strcmp(vars[i].name, "isothermal")==0)      isothermal = (int) vars[i].value;
        if (strcmp(vars[i].name, "unsaturated")==0)     unsaturated = (int) vars[i].value;
        if (strcmp(vars[i].name, "torque_to_acc")==0)   torque_to_acc = (int) vars[i].value;
        if (strcmp(vars[i].name, "read_sig_from_file")==0)   read_sig_from_file = (int) vars[i].value;
        if (strcmp(vars[i].name, "stellar_torque")==0)  stellar_torque = (int) vars[i].value;
        if (strcmp(vars[i].name, "R_star_to_sun")==0)   R_star_to_sun = vars[i].value * 0.00465047; //in au
        if (strcmp(vars[i].name, "R_co")==0)            R_co = vars[i].value;
    }    
    
    // Check if the disc is really meant to be isothermal
    if (isothermal == 0 && ADIABATICINDEX != 1) 
        reb_exit("Your disc supposed to be locally-isothermal but the adiabatic index is not 1.\n");
    
  
}


int find_planets_number(char *which_input){
    int counter=0;
    // planet 0 is the primary/star
    if (strcmp(which_input,"fargo") == 0){
        DIR *dir_here = opendir(".");
        struct dirent *dir;
        int norb = 0, nbpl = 0;
        if (dir_here){
            while ((dir=readdir(dir_here)) != NULL){
                if (*(dir->d_name) != '.' ){
                    char *isorbit = strstr(dir->d_name, "orbit");
                    if (isorbit != NULL){  
                        norb++;
                    }
                    char *isbpl = strstr(dir->d_name, "bigplanet");
                    if (isbpl != NULL){
                        nbpl++;
                    }
                }
            }
            if (norb != nbpl){
                printf("norb=%d, npl=%d\n", norb, nbpl);
                reb_exit("Numebr of orbit.dat and bigplanet.dat files are not equal. Some files are missing.\n");
            } else
                counter = norb;
        }
        closedir(dir_here);
        tstart = 0;
    } else if ((strcmp(which_input,"rebound") == 0) || (strcmp(which_input,"restartdisc") == 0)){
        DIR *dir_here = opendir(".");
        struct dirent *dir;
        if (dir_here){
            while ((dir=readdir(dir_here)) != NULL){
                char *iselement = strstr(dir->d_name, "restart_final.bin");
                if (iselement != NULL){
                    struct reb_simulation* r = reb_create_simulation_from_binary("restart_final.bin");
                    counter = r->N-1;
                    tstart = r->t;
                    reb_free_simulation(r);
                }
            }
            if (counter == 0)
                reb_exit("restart_final.bin does not exist.\n");
        }
        closedir(dir_here);
    } else if (strcmp(which_input,"disc") == 0) { 
        FILE *input;
        char s[512], filename[512];
        // opening the planets.par file. It must include 
        // planet's name, mass, semi-major axis, eccentrcity, migration type, tau_a, tau_e
        // in a single line per planet
        sprintf (filename, "planets.par");
        input = fopen (filename, "r");
        if (input == NULL){
            printf("I need a file named 'planets.par'. It must include 6 columns\n \
                    planet's name, Mass, Semi-majorAxis, Eccentricity, Migration type, tau_a, tau_e,\
                     in a single line per planet. \n");
            reb_exit("I can't find planets.par file :-(.\n");
        }
        while (fgets(s, 510, input) != NULL) {
            if (isalnum(s[0])) // check if first character alphabetical or numerical
                counter++;
        }
        fclose (input);
        tstart = 0;
    } else {
        printf("Invalid input paramter in find_planets_number...\n");
        exit(1);
    }
    
    printf("Total number of planets = %d\n", counter);
    return counter;
}

void read_planets(init_planet *planets, char *which_input, int nplanets){
    float mass, a, e, taua, taue, tpert, apert;
    if (strcmp(which_input,"fargo") == 0){
        float foo, paromega, xpl, ypl, theta;
        int ifoo;
        FILE *input;
        char filename[50], line[610];
        for (int counter=0; counter<nplanets; counter++){
            // Reading the mass and position of the planets
            sprintf(filename, "bigplanet%d.dat",counter);
            input = fopen(filename, "r");
            while(!feof(input))
                fgets(line, 600, input);
            sscanf(line, "%d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %d %e %e %e %e", \
                &ifoo, &xpl, &ypl, &foo, &foo, &mass, &foo, &foo, &foo, &foo, &foo, &foo, &foo, &foo, &foo, &ifoo, \
                &foo, &foo, &foo, &foo);
            fclose(input);
            // Reading the orbital elements
            sprintf(filename, "orbit%d.dat",counter);
            input = fopen(filename, "r");
            while(!feof(input))
                fgets(line, 600, input);
            sscanf(line, "%e %e %e %e %e %e %e %e", &foo, &e, &a, &foo, &foo, &paromega, &foo, &foo);
            theta = atan2(ypl, xpl);
            planets[counter].a = a;
            planets[counter].mass = mass;
            planets[counter].e = e;
            strcpy(planets[counter].mig_type, "0");
            planets[counter].tau_a = 0;
            planets[counter].tau_e = 0;
            planets[counter].inc = 0;
            planets[counter].Omega = 0;
            planets[counter].omega = paromega;
            planets[counter].f = theta-paromega;
            printf("Planet %d: semi-major-axis = %lg, mass = %lg, eccentricity = %lg \n", \
                           counter, planets[counter].a, planets[counter].mass, planets[counter].e);
            fclose(input);
        }

    } else if ((strcmp(which_input,"rebound") == 0) || (strcmp(which_input,"restartdisc") == 0)){
        printf("Restarting from the final rebound archive, if you want some time else, you need to modify the filename\n");
        struct reb_simulation* r = reb_create_simulation_from_binary("restart_final.bin");
        if (r == NULL)
            reb_exit("Problem in openning the archive file.\n");
        struct reb_particle *particles = r->particles;
        for (int counter=1; counter < r->N; counter++){
            struct reb_orbit o = reb_tools_particle_to_orbit(r->G, particles[counter], particles[0]);
            planets[counter-1].a = o.a;
            planets[counter-1].mass = particles[counter].m;
            planets[counter-1].e = o.e;
            strcpy(planets[counter-1].mig_type, "0");
            planets[counter-1].tau_a = 0;
            planets[counter-1].tau_e = 0;
            planets[counter-1].inc = o.inc;
            planets[counter-1].Omega = o.Omega;
            planets[counter-1].omega = o.omega;
            planets[counter-1].f = o.f;
            if (strcmp(which_input,"restartdisc") == 0) {
                char line[352], name[100], migtype[20];    
                int counter=0;   
                FILE *input;
                input = fopen("planets.par", "r");
                if (input == NULL)
                    reb_exit("Problem in opening planets.par file :-(.\n");
                while( fgets(line, 350, input) != NULL){
                    if (isalnum(line[0])){
                        sscanf(line, "%s %e %e %e %s %e %e %e %e", name, &mass, &a, &e, migtype, &taua, &taue, &tpert, &apert);
                        strcpy(planets[counter].mig_type, migtype);
                        planets[counter].tau_a = taua;
                        planets[counter].tau_e = taue;
                        planets[counter].t_pert = tpert;
                        planets[counter].a_pert = apert;
                        counter++;
                    }
                }
                fclose (input);
            }
            printf("Planet %d: semi-major-axis = %lg, mass = %lg, eccentricity = %lg, migration_type=%s \n", \
                           counter-1, planets[counter-1].a, planets[counter-1].mass, planets[counter-1].e, planets[counter-1].mig_type);
        }
        reb_free_simulation(r);
    } else if (strcmp(which_input,"disc") == 0) { 
        char line[352], name[100], migtype[20];    
        int counter=0;   
        FILE *input;
        input = fopen("planets.par", "r");
        if (input == NULL)
            reb_exit("Problem in opening planets.par file :-(.\n");
        while( fgets(line, 350, input) != NULL){
            if (isalnum(line[0])){
                sscanf(line, "%s %e %e %e %s %e %e %e %e", name, &mass, &a, &e, migtype, &taua, &taue, &tpert, &apert);
                planets[counter].a = a;
                planets[counter].mass = mass;
                planets[counter].e = e;
                strcpy(planets[counter].mig_type, migtype);
                planets[counter].tau_a = taua;
                planets[counter].tau_e = taue;
                planets[counter].inc = 0;
                planets[counter].Omega = 0;
                planets[counter].omega = 0;
                planets[counter].f = 0;
                planets[counter].t_pert = tpert;
                planets[counter].a_pert = apert;
                printf("Planet %d: semi-major-axis = %lg, mass = %lg, eccentricity = %lg, migration type=%s \n", \
                               counter, planets[counter].a, planets[counter].mass, planets[counter].e, planets[counter].mig_type);
                counter++;
            }
        }
        fclose (input);
    }   
}


void additional_forces(struct reb_simulation* const r){
    double r2, dist, h, sigma_p, beta_temperature, omega, e, a, e_hat=0;
    int i, n_objects;
    double Q, ax=0, ay=0, vx, vy, x, y, mass, ksi, K;
    double gamq, gamq2, denominator, Gammaeff, xs, Pk, Pnu;

    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0]; // calculate migration forces with respect to center of mass;
    
    //sg_wrong_factor is the factor that makes the same torque as in hydro if SG is ignored

    //Number of objects
    n_objects = r->N;
    
    for (i=1; i<n_objects; i++){
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, particles[i], particles[0]);
        e = o.e;
        a = o.a;

        x = particles[i].x - com.x;
        y = particles[i].y- com.y;
        vx = particles[i].vx- com.vx;
        vy = particles[i].vy- com.vy;
        mass = particles[i].m;
        
        double G = r->G;
        double mu = G*(com.m + mass);

        r2 = x*x+y*y;
        
        dist = sqrt(r2);
        if (dist < 0.05) //0.05 is the inner boundary
            reb_exit("One of the particles gets very close to the star.\n");
        
        h = calculate_aspect_ratio(a); //(isothermal) aspect ratio = cs_iso/vk
        sigma_p = calculate_sigma_at_rp(a); 
        //sigma_p *= -0.062*(-19.45+pow(time_now,0.26)); //decay for the surface density profile
        alpha_sigma = calculate_sigma_slope(a); 
        beta_temperature = calculate_temperature_slope(a);
        omega = sqrt(mu/a/a/a);
        
        // torque scaling coeficients, see section 3 and 5.6 of the paper
        // These torque components are only for locally isothermal models.
        Gamma0[i] = mass/h/h * a * sigma_p;
        
        
        if (isothermal != 0){
	    reb_exit("Under construction :-(.\n");
            Q = 2*DIFFUSIVITY/3/(h*h*h)/sqrt(dist);
            gamq = ADIABATICINDEX*Q;
            gamq2 = gamq*gamq;
            denominator = sqrt(pow(gamq2+1,2)-16*(Q*Q)*sqrt(ADIABATICINDEX-1)) + gamq2 - 1;
            denominator *= 2;
            denominator = 0.5*sqrt(denominator);
            denominator += gamq;
            Gammaeff = 2*gamq/denominator;
            xs = 1.1/pow(Gammaeff,0.25)*pow(0.4/THICKNESSSMOOTHING,0.25)*sqrt(mass/h);
            Pk = sqrt(sqrt(dist)*xs*xs*xs/2/M_PI/DIFFUSIVITY);
        } else {  //if locally isothermal
            Gammaeff = ADIABATICINDEX;
            Pk = 0;
        }
        
        double visc = caclulate_viscosity(dist, h);
        //Here I added the cavity viscosity reduction, it is only used when the torque is partially-saturated
        double cavity_radius = 1.0, cavity_width = 0.6, cavity_ratio = 100.0;
        double rmin = cavity_radius - cavity_width*h;
        double rmax = cavity_radius + cavity_width*h;
        if ((dist >= rmin) && (dist <= rmax))
	    visc *= exp((rmax-dist)/(rmax-rmin)*log(cavity_ratio));
        if (dist<rmin)
            visc *= cavity_ratio;
        K = sqrt(dist) /(2*M_PI*visc) ;
        xs = 1.1/pow(Gammaeff,0.25)*pow(0.4/THICKNESSSMOOTHING,0.25)*sqrt(mass/h);
        ksi = beta_temperature - (Gammaeff-1)*alpha_sigma;
        Pnu = 2./3. * sqrt(xs*xs*xs*K);
        
        // Barotropic vortensity related horseshoe drag Eq. 4 with modification from 2010 paper
        GammaHSB[i] = 1.1 * (1.5-alpha_sigma) * (0.4/THICKNESSSMOOTHING);
        // Entropy related linear corotation torque Eq. 7 with modification from 2010 paper
        GammaLE[i] = 2.2 * ksi * pow(0.4/THICKNESSSMOOTHING,0.71) - 1.4 * ksi/Gammaeff * pow(0.4/THICKNESSSMOOTHING,1.26);
        // Barotropic vortensity related linear corotaion torque. it is zero in LISO models
        GammaLB [i]= 0.7 * (1.5-alpha_sigma) * pow(0.4/THICKNESSSMOOTHING,1.26);
        //Entropy related horseshoe drag, Eq. 5 with modifications from Paardekooper 2010 paper. This is also zero for LISO
        GammaHSE[i] = ksi/Gammaeff *(0.4/THICKNESSSMOOTHING) * (10.1*sqrt(0.4/THICKNESSSMOOTHING)-2.2);
        // Lindblad torque relation 47 
        GammaL[i] = (-2.5 - 1.7*beta_temperature + 0.1*alpha_sigma)* pow(0.4/THICKNESSSMOOTHING,0.71);
        
        if (unsaturated == 0){
            if (isothermal != 0)
                GammaCTotal[i] = GammaHSE[i] + GammaLE[i] + GammaHSB[i] + GammaLB[i];
            else  //if locally isothermal 
                GammaCTotal[i] = GammaLE[i] + GammaHSB[i];        
        } else {
            GammaCTotal[i] = GammaHSE[i] * Ffunc(Pnu)*Ffunc(Pk) * sqrt(Gfunc(Pnu)*Gfunc(Pk));
            GammaCTotal[i] += GammaLE[i] * sqrt((1-Kfunc(Pnu)) * (1-Kfunc(Pk)));
            GammaCTotal[i] += GammaHSB[i] * Gfunc(Pnu) * Ffunc(Pnu);
            GammaCTotal[i] += GammaLB[i] * (1-Kfunc(Pnu));
        }
        
        
        // Applying the eccentricity correction on the corotation torque
        // using the relation 8 and 10 of Fendyke 2013 
        e_hat = e/h;
        if (without_e_damping == 0){
            GammaCTotal[i] *= exp(-1 * e/ (0.5*h+0.01)); 
        }
        
        if (without_e_damping == 2){
            GammaCTotal[i] *= exp(-1 * e/ (0.5*h+0.01)); 
            double p_e = (1+sqrt(0.444*e_hat)+pow(0.352*e_hat,6))/(1-pow(0.495*e_hat,4));
            GammaL[i] /= p_e;
	}
        if (without_e_damping == 3){
            GammaL[i] /= (1-pow(0.495*e_hat,4))/(1+sqrt(0.444*e_hat)+pow(0.352*e_hat,6));
        }
        
        // The below term is added to simulate the torque for a moving planet in hydro simulations w/o Axi-SG term.
        if (without_sg_term == 1){
            Q = h/(M_PI*dist*dist*sigma_p);
            double log_Qh = log10(Q*h);
            if (alpha_sigma == 1.5 && flaring == 0.5){ 
                //sg_wrong_factor = (1+0.6 * pow(Q*h,-1.02)) * 1.034; //This last coefficient is for the correction from PBK to the fixed models
                double wrong_fact = pow(10,0.13*log_Qh*log_Qh - 0.58*log_Qh -0.64);
                sg_wrong_factor[i] = 1+ wrong_fact;
            } else if (alpha_sigma == 0.5 && flaring == 0.0){
                //sg_wrong_factor = (1+0.1 * pow(Q*h,-0.64)) * 1.06;
                double wrong_fact = pow(10,0.22*log_Qh*log_Qh - 1.09*log_Qh -0.17);
                sg_wrong_factor[i] = 1+ wrong_fact;
            } else {
                reb_exit("I don't have the correction term for this model. \n");
            }
            GammaL[i] *= sg_wrong_factor[i];
        }
        
        Gamma[i] = GammaCTotal[i] + GammaL[i];
        Gamma[i] *= Gamma0[i]/Gammaeff;
                
        // Converting the torque to acceleration
        //positive eccentricity damping = damping
        int i_p = i-1;
        if (strcmp(planets[i_p].mig_type,"torque") == 0) {
            if (without_e_damping == 2)
                tau_ecc[i] =  1.282*(1-0.14*e_hat*e_hat+0.06*pow(e_hat,3)) *h*h*a*a*omega / Gamma0[i]; //From Cresswell+2008
            else
                tau_ecc[i] =  2.6 * h*h * dist*dist * omega / Gamma0[i]; //From Baruteau+2014
            //tau_mig positive = inward migration        
    	    tau_mig[i] = - sqrt(o.a*(1-e*e)) / Gamma[i]; 
            tau_a[i] = 0.5/(1./tau_mig[i]+e*e/(1-e*e)/tau_ecc[i]);
        } else if (strcmp(planets[i_p].mig_type,"timescale") == 0) {
            tau_a[i] =  planets[i_p].tau_a* 2*M_PI;
            tau_ecc[i] = planets[i_p].tau_e * 2*M_PI;
            double inverse = 1./2/(tau_a[i])-e*e/(1-e*e)/tau_ecc[i];
            tau_mig[i] = 1./inverse;
        } else if (strcmp(planets[i_p].mig_type,"oneside") == 0) {
            // Hand-made planetary trap one sided
            tau_a[i] =  planets[i_p].tau_a* 2*M_PI;
            tau_ecc[i] = planets[i_p].tau_e * 2*M_PI;
            double r_o = 1.15, delta = 0.05;
            double f_corr = tanh((a-r_o)/delta);
            tau_a[i] /= f_corr;
            double inverse = 1./2/(tau_a[i])-e*e/(1-e*e)/tau_ecc[i];
            tau_mig[i] = 1./inverse;
        } else if (strcmp(planets[i_p].mig_type,"twoside") == 0) {
            // Hand-made planetary trap two sided
            tau_a[i] =  planets[i_p].tau_a* 2*M_PI;
            tau_ecc[i] = planets[i_p].tau_e * 2*M_PI;
            double r_o = 1.15, delta = 0.05, r_i=0.93, depth=0.9;
            //double f_corr = 1+depth*(tanh((a-r_o)/delta)+tanh((r_i-a)/delta));
            double scale = (1-depth);
            double f_corr = 1.5+(tanh((a-r_o)/delta)+depth*tanh((r_i-a)/delta))/scale;
            tau_a[i] /= f_corr;
            double inverse = 1./2/(tau_a[i])-e*e/(1-e*e)/tau_ecc[i];
            tau_mig[i] = 1./inverse;
        } else if (strcmp(planets[i_p].mig_type,"withbump") == 0) {
            // Hand-made double planetary trap two-sided with bump
            tau_a[i] =  planets[i_p].tau_a* 2*M_PI;
            tau_ecc[i] = planets[i_p].tau_e * 2*M_PI;
            double delta = 0.2, depth=2., r_ref=0.49, frac=0.4, shift=0.8;
            double f_corr = frac*tanh((a-r_ref)/delta);
            f_corr += depth*((a-(1+shift*delta)*r_ref)/delta)*exp(-pow((a-r_ref)/delta,2));
            f_corr += 1-frac;
            tau_a[i] /= f_corr;
            double inverse = 1./2/(tau_a[i])-e*e/(1-e*e)/tau_ecc[i];
            tau_mig[i] = 1./inverse;
        } else {
            reb_exit("Not a valid migration type\n");
        }

        if (stellar_torque == 0){
            GammaStar[i] = 1.5 * mass * pow(R_star_to_sun,5) * pow(a,-6) * pow(10,-3.25);
            if (a < R_co)
                GammaStar[i] *= -1;
            double gamma_new = - sqrt(o.a*(1-e*e)) / tau_mig[i] + GammaStar[i]; 
            tau_mig[i] = - sqrt(o.a*(1-e*e)) / gamma_new; 
            tau_a[i] = 0.5/(1./tau_mig[i]+e*e/(1-e*e)/tau_ecc[i]);
        }

        if (torque_to_acc == 0){
            // The method of Rebound
            double hh = x*vy - y*vx;
            //double v = sqrt ( vx*vx + vy*vy);
            double r = dist;
            //double vr = (x*vx + y*vy)/r;
            //double ex = 1./mu*( (v*v-mu/r)*x - r*vr*vx );
            //double ey = 1./mu*( (v*v-mu/r)*y - r*vr*vy );
            //double e = sqrt( ex*ex + ey*ey);        // eccentricity
            //double a = -mu/( v*v - 2.*mu/r );            // semi major axis
            double prefac1 = 1./(1.-e*e) /tau_ecc[i]/1.5;
            double prefac2 = 1./(r*hh) * sqrt(mu/a/(1.-e*e))  /tau_ecc[i]/1.5;
            ax = -1. * vx/tau_mig[i] - vx*prefac1 - (hh*y)*prefac2;
            ay = -1. * vy/tau_mig[i] - vy*prefac1 +(hh*x)*prefac2;
        }
        if (torque_to_acc == 2){
            // The method of Papaloizou
            double vdotr = vx*x + vy*y;
            ax =  -1. * vx/tau_mig[i] -2./dist/dist * vdotr * x/tau_ecc[i];
            ay =  -1. * vy/tau_mig[i] -2./dist/dist * vdotr * y/tau_ecc[i];
        }
        
        // Updating the acceleration 
        particles[i].ax += ax;
        particles[i].ay += ay;
        com = reb_get_com_of_pair(com,particles[i]);
        
    }
}


double Ffunc(double p){
    return 1./(1+p*p/1.3/1.3);
}

double Gfunc(double p){
    double crit;
    crit = 8./45/M_PI;
    if (p < sqrt(crit))
        return 16./25.*pow(crit,-3./4)*pow(p,1.5);
    else 
        return 1-9./25.*pow(crit,4./3)*pow(p,-8./3);
}

double Kfunc(double p){
    double crit;
    crit = 28./45./M_PI;
    if (p < sqrt(crit))
        return 16./25.*pow(crit,-3./4)*pow(p,1.5);
    else 
        return 1-9./25.*pow(crit,4./3)*pow(p,-8./3);
}


void heartbeat(struct reb_simulation* const r){
    // Output some information to the screen and into a binary archive file 
    if(reb_output_check(r, tdump)){
        reb_output_timing(r, tmax);
        char name[50];
        FILE *filename;
        struct reb_particle* particles = r->particles;
        for (int i=1; i<r->N; i++){
            //print orbital elements of planets separately
            sprintf(name, "orbit%d.txt", i);
            filename = fopen(name,"a");
            struct reb_orbit o = reb_tools_particle_to_orbit(r->G, particles[i], particles[0]);
            fprintf(filename,"%e\t%e\t%e\t%e\t%e\n",tstart+r->t,o.e, o.a, o.M, o.pomega); //time, semi-major axis, eccentricity
            fclose(filename);
            //print torques separately
            sprintf(name, "torque%d.txt", i);
            filename = fopen(name,"a");
            fprintf(filename, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", tstart+r->t, Gamma[i], Gamma0[i], \
                GammaL[i], GammaCTotal[i], GammaHSB[i], GammaHSE[i], GammaLE[i], GammaLB[i], sg_wrong_factor[i], GammaStar[i]);
            fclose(filename);
            //print timescales of planets separately
            sprintf(name, "tau%d.txt", i);
            filename = fopen(name,"a");
            fprintf(filename, "%e\t%e\t%e\t%e\n", tstart+r->t, tau_mig[i], tau_a[i], tau_ecc[i]);
            fclose(filename);
            if(planets[i-1].t_pert > 1e-10){
                if(reb_output_check(r, planets[i-1].t_pert)){
                    //Make some perturbation on the velocity of the planets
                    double vx = particles[i].vx, vy = particles[i].vy;
                    double random = (2*(rand()/(double)RAND_MAX) - 1)*planets[i-1].a_pert;
                    //printf("%d, %e, %e\n", i, random, particles[i].m);
                    particles[i].vx += vx* random;
                    particles[i].vy += vy* random;
                }
            }
        }
    }
    //Make 10 binary simulation archives per simulations
    if(reb_output_check(r, tmax/10.)){
        char filename[50];
        sprintf(filename, "restart%d.bin", counter_archive);
        reb_output_binary(r, filename);
        counter_archive++;
    }
}

double calculate_aspect_ratio(double dist){
    //If the disc is locally isothermal, the aspect ratio is calculated simply by h0 and flaring index
    //otherwise a function must be added for calculating the aspect ratio in a radiative disc
    return  h0*pow(dist,flaring);
}

double calculate_sigma_at_rp(double dist){
    // Calculating the surface density at the location of the planet
    // needs to be modified for a non-powerlow profile
    int i=0;
    if (read_sig_from_file != 0){
        return sigma0*pow(dist, -alpha_sigma);
    } else {
        while(rad_array[i] < dist)
            i++;
        return sigma_array[i];
    }
}

double calculate_sigma_slope(double dist){
    // Calculating the surface density slope (dlnsigma/dlnr) at the location of the planet
    // needs to be modified for a non-powerlow profile
    int i=0, tolerance=1;
    if (read_sig_from_file != 0){
        return alpha_sigma;
    } else {
        while (rad_array[i] < dist)
            i ++;
        double slope = log(sigma_array[i+tolerance])-log(sigma_array[i-tolerance]);
        slope /= log(rad_array[i+tolerance])-log(rad_array[i-tolerance]);
        /*printf("r_p=%lg, sigma_p=%lg, r-1=%lg, sigma-1=%lg, r1=%lg, sigma1=%lg, slope=%f\n", rad_array[i], sigma_array[i], \
        rad_array[i-tolerance], sigma_array[i-tolerance], rad_array[i+tolerance], sigma_array[i+tolerance], slope);*/
        return -slope;       
    }
} 

double calculate_temperature_slope(double dist){
    // Calculating the temperature density slope (dlnT/dlnr) at the location of the planet
    // needs to be modified for a non-powerlow profile
    return -2*flaring+1; 
}

double caclulate_viscosity(double dist, double h_iso){
    // Calculating the viscosity the location of the planet
    return alpha_visc*h_iso*h_iso*sqrt(ADIABATICINDEX)*sqrt(dist);
}

void read_sigma(){
    char filename[512], line[352];
    int counter=0;   
    FILE *input;
    
    sprintf(filename, "sigma.dat");
    input = fopen(filename, "r");
    if (input == NULL){
        reb_exit("Problem in opening sigma.dat file :-(.\n \
                  You need a file named sigma.dat with the first column\n\
                  being the radii and the second surface density.");
    }
    while( fgets(line, 350, input) != NULL){
        if (isdigit(line[0])){
            sscanf(line, "%lf  %lf", &rad_array[counter], &sigma_array[counter]);
            counter++;
        }
    }
    printf("\nYou have %d points in the radial direction.\n", counter);
    fclose (input);
}

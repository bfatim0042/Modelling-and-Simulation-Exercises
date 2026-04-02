////#-- ========================================================================================
//#--                Part 2, Exercise 1 - Bea Fatima, Student ID: 1725181
//#-- ========================================================================================
//THIS CODE CALCULATES PE AND KE FOR A LENNARD-JONES SYSTEM USING THE VERLET ALGORITHM AND ANDERSEN THERMOSTAT

//Using Exercise 8 with LJ interactions as our base

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#define NDIM 3
#define N 512

/* Initialization variables */
const int    nsteps = 10000; //time steps
const double density = 0.5; //density
const double r_cut   = 2.5; //cutoff distance
const double dt      = 0.001; //time step
const double init_temp   = 1.5; //initial temp
const double nu = 1.0;   // collision frequency
const char*  init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0; 
double box[NDIM]; //box dim
double e_cut; //potential energy cutoff

double x[N][NDIM];   // positions
double v[N][NDIM];   // velocities
double f[N][NDIM];   // forces
double xm[N][NDIM];  // previous positions (used for Verlet)


void init(){ //adapted from pseudocode (Table 1.1, pg 2)

    double sumv[NDIM] = {0}; //total velocity in each dim
    double sumv2 = 0.0; //total velocity squared

    for(int i=0;i<n_particles;i++){ //loop over particles
        for(int d=0;d<NDIM;d++){ //loop over dimensions
            v[i][d] = 2.0*(dsfmt_genrand() - 0.5); //generate random velocity between -1 and 1
            sumv[d] += v[i][d]; //sum velocities for total momentum
            sumv2   += v[i][d]*v[i][d]; //sum velocities squared
        }
    }

    for(int d=0;d<NDIM;d++) //average velocity in each dim
        sumv[d] /= n_particles;

    sumv2 /= n_particles; //average velocity squared

    double fs = sqrt(3.0*init_temp / sumv2); //scaling factor

    for(int i=0;i<n_particles;i++){ //loop over particles 
        for(int d=0;d<NDIM;d++){//loop over dimensions
            v[i][d] = (v[i][d] - sumv[d]) * fs; //add scaling factor
            xm[i][d] = x[i][d] - v[i][d]*dt; //updates previous positions for Verlet
        }
    }
}

void force(double *en){ //adapted from pseudocode (Table 1.2, pg 3)

    *en = 0.0; //set total energy to 0

    for(int i=0;i<n_particles;i++) //loop over particles
        for(int d=0;d<NDIM;d++) //loop over dimensions
            f[i][d] = 0.0; //set forces to 0

    for(int i=0;i<n_particles-1;i++){ //loop over particle pairs
        for(int j=i+1;j<n_particles;j++){

            double xr[NDIM];
            double r2 = 0.0;

            for(int d=0;d<NDIM;d++){ //loop over dim

                xr[d] = x[i][d] - x[j][d]; //distance between particles in each dim

                //periodic boundary conditions
                if(xr[d] >  0.5*box[d]) xr[d] -= box[d]; 
                if(xr[d] < -0.5*box[d]) xr[d] += box[d];

                r2 += xr[d]*xr[d]; //distance square
            }

            if(r2 < r_cut*r_cut && r2 > 1e-12){ //if particles are within cutoff distance

                double r2i = 1.0/r2; //(1/dist^2)
                double r6i = r2i*r2i*r2i; //(1/dist^6)
                double ff = 48.0*r2i*r6i*(r6i - 0.5); //force factor

                for(int d=0;d<NDIM;d++){
                    f[i][d] += ff*xr[d]; //update force for particle i
                    f[j][d] -= ff*xr[d]; //update force for particle j
                }

                *en += 4.0*r6i*(r6i - 1.0) - e_cut; //update potential total energy
            }
        }
    }
}


double gauss(double sigma){  //generating Gaussian random number (box-muller method)
    double u1 = dsfmt_genrand();
    double u2 = dsfmt_genrand();
    return sigma * sqrt(-2.0*log(u1)) * cos(2.0*M_PI*u2); 
}

double compute_energy(double en){ //adapted from pseudocode (Table 1.3, pg 5)
    double sumv2 = 0.0;
    for(int i = 0; i < n_particles; i++){
        for(int d = 0; d < NDIM; d++){
            sumv2 += v[i][d]*v[i][d]; 
        }
    }
    return (en + 0.5*sumv2)/n_particles; //total energy per particle
}


//referencing Table 1.5, pg 7 for Andersen thermostat pseudocode structure

void integrate(int switch_step, double *en, double *temp){

    if(switch_step == 1){ //first step

        for(int i = 0; i < n_particles; i++){ //loop over particles
            for(int d = 0; d < NDIM; d++){ //loop over dimensions
                x[i][d] = x[i][d] + dt*v[i][d] + 0.5*dt*dt*f[i][d]; //verlet position update
                v[i][d] = v[i][d] + 0.5*dt*f[i][d]; //verlet velocity update
            }
        }

    } else if(switch_step == 2){ //second step

        double tempa = 0.0; 
        for(int i = 0; i < n_particles; i++){ //loop over particles
            for(int d = 0; d < NDIM; d++){ //loop over dimensions
                v[i][d] = v[i][d] + 0.5*dt*f[i][d]; //update velocity
                tempa += v[i][d]*v[i][d]; //sum of velocity squared to calculate temperature
            }
        }
        *temp = tempa / (3.0 * n_particles); //instantaneous temperature
        double sigma = sqrt(init_temp); //andersen heat bath

        for(int i = 0; i < n_particles; i++){

            if(dsfmt_genrand() < nu*dt){ //test for collision with heat bath
                for(int d = 0; d < NDIM; d++){
                    v[i][d] = gauss(sigma); //give particle Gaussian velocity
                }
            }
        }
    }
}

void read_data(void){ //same as in exercise 8
    FILE* fp = fopen(init_filename, "r");
    int n, d;
    double dmin, dmax;
    fscanf(fp, "%d", &n_particles);
    for(d = 0; d < NDIM; ++d){
        fscanf(fp, "%lf %lf", &dmin, &dmax);
        box[d] = fabs(dmax - dmin);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fscanf(fp, "%lf\t", &x[n][d]);
        double diameter;
        fscanf(fp, "%lf\n", &diameter);
    }
    fclose(fp);
}

void set_density(void){ //same as in exercise 8
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = n_particles / density;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) x[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}

int main(){

    size_t seed = time(NULL);
    dsfmt_seed(seed);

    read_data(); //read initial config from fcc.dat
    set_density(); //rescale density

    double r6_cut = pow(1.0 / r_cut, 6.0); //cutoff distance
    e_cut = 4.0 * (r6_cut*r6_cut - r6_cut); //potential energy cutoff

    for(int d = 0; d < NDIM; ++d)
        assert(r_cut <= 0.5 * box[d]); //make sure cutoff is less than half the box size

    init();

    double en; 
    force(&en);   // initial force

    //from initial code
    double volume = 1.0;
    for(int d = 0; d < NDIM; ++d)
        volume *= box[d]; 

    printf("Starting volume: %f\n", volume);
    printf("Starting seed: %lu\n", seed);

    FILE* fp = fopen("energy.dat", "w");

    double t = 0.0;
    double temp;

    for(int step = 0; step < nsteps; ++step){ //main simulation loop

        integrate(1, &en, &temp);   //first part of e.o.m
        force(&en);                 //update forces and potential energy
        integrate(2, &en, &temp);   //second part of e.o.m

        double sumv2 = 0.0;
        for(int i = 0; i < n_particles; i++){
            for(int d = 0; d < NDIM; d++){
                sumv2 += v[i][d]*v[i][d]; //sum of velocity squared for kinetic energy
            }
        }

        double ke = 0.5 * sumv2 / n_particles; //kinetic energy per particle
        double pe = en / n_particles; //potential energy per particle

        fprintf(fp, "%f\t%f\t%f\n", t, pe, ke); //write time, potential energy, kinetic energy

        t += dt; //update time
    }

    fclose(fp);
    return 0;
}

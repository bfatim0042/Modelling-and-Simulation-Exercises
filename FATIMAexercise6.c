//#-- ========================================================================================
//#--                Exercise 6 - Bea Fatima, Student ID: 1725181
//#-- ========================================================================================

#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 10000;
const int output_steps = 100;
const double packing_fraction = 0.05;
const double diameter = 1.0;
const double delta = 0.1;
const char* init_filename = "liquidfcc.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];


#define MAXBIN 200
double hist[MAXBIN]; //histogram to store particle pair dist count
double gr[MAXBIN]; // distribution function g(r)
double delr; //width of bin 

/* Functions */
void read_data(void){
    /*--------- Your code goes here -----------*/
//open file and read data, error message if it is null
    FILE *fccdata = fopen(init_filename, "r");
    if (fccdata == NULL){
        printf("Could not open file, is null");
        return;
    }

//check number of particles and store into n_particles
fscanf(fccdata, "%d", &n_particles);

//check dimensions and store into box array
//only storing max since min = 0
for (int d = 0; d < NDIM; d++) {
    double min, max;
    fscanf(fccdata, "%lf %lf", &min, &max);
    box[d] = max;

}
//read coordinates and store into array
//loop over all particles
for (int n = 0; n < n_particles; ++n) {
    //loop over all dimensions
    for (int d = 0; d < NDIM; ++d) {
        //read double and store into r[n][d] 2-d array (n being particle # and d being which dim x,y,or z)
        fscanf(fccdata, "%lf", &r[n][d]);
    }
    //read and discard diameter value so it doesnt count
    double skip;
    fscanf(fccdata, "%lf", &skip);
}
//close file
fclose(fccdata);

//verifying data amounts
//checking # of particles
printf("n_particles = %d\n", n_particles);

//checking box dimensions
printf("Box: %lf %lf %lf\n", box[0], box[1], box[2]);

//coords for first particle
printf("\nFirst particle: ");
for (int d = 0; d < NDIM; ++d) printf("%lf ", r[0][d]);
    printf("\n");
}

//code from the previous exercise
int move_particle(void){
//pick random particle index out of n_particles
int n = (int)( dsfmt_genrand() * n_particles );

//random displacement between -delta to delta in each direction
//uses random number generator to get number between [-1,1], then multiply by delta
double newpos[NDIM];
//loop over dimensions so we can add displacement to x,y,z
for(int d = 0; d < NDIM; d++){
    double displacement = (2.0 * dsfmt_genrand() - 1.0) * delta;
    newpos[d] = r[n][d] + displacement;
   
}
 //make sure particles stay in box, periodic boundary cond.
 //loop over dimensions to see if it is out of bounds
 for(int d = 0; d < NDIM; d++){
    //if out of bounds, fix by wrapping it around by length of box
        //if coord less than 0, add box length 
        if(newpos[d] < 0.0)
            newpos[d] += box[d];
        //if coord greater than 0, subtract box length 
        else if(newpos[d] >= box[d])
            newpos[d] -= box[d];
    }

//check for overlapping, and reject in that case
for (int m = 0; m < n_particles; m++) {

    //case where it checks with itself
    if (m == n) 
    continue;

    double dist = 0.0;
    for (int d = 0; d < NDIM; d++) {
        double dx = newpos[d] - r[m][d];
//use nearest image convention 
//distances between particles in any of the three directions are always larger than −L/2 and always smaller than L/2.
            if(dx >  0.5 * box[d]) dx -= box[d];
            if(dx < -0.5 * box[d]) dx += box[d];
            dist += dx * dx; //squared distance
        }
        //if squared distance is less than squared diameter, it overlaps
        //reject
        if(dist < (diameter * diameter)){
            return 0;   //reject move
        }
    }
    //no overlaps, accept the newpos with displacement as actual
    for(int d = 0; d < NDIM; d++){
        r[n][d] = newpos[d];
    }
    //change is accepted
    return 1;
}
void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords_step%07d.dat", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%f\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

void set_packing_fraction(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
}
//taken from pseudocode in Lecture 4 slides
void gofr(){

    int i, j, d;

    // loop over particle pairs
    for(i = 0; i < n_particles-1; i++){

        for(j = i+1; j < n_particles; j++){

            double dist = 0.0;
// find squared distance between particles i and j
            for(d = 0; d < NDIM; d++){

                double dx = r[i][d] - r[j][d];

                // apply periodic boundary conditions
                // nearest image convention
                if(dx >  0.5 * box[d]) dx -= box[d];
                if(dx < -0.5 * box[d]) dx += box[d];

                dist += dx * dx;
            }
// squared distance -> actual distance
            double rij = sqrt(dist);
// find correct histogram bin 
            int bin = (int)(rij / delr);
 // add to histogram if within range
            if(bin < MAXBIN){
                hist[bin] += 2.0;
            }
        }
    }
}
void normalise(){

     // average histogram over all steps
    for(int bin = 0; bin < MAXBIN; bin++){
        hist[bin] /= mc_steps;
    }

    // calculate number density 
    double volume = box[0]*box[1]*box[2];
    double rho = n_particles / volume;

    double constant = 4.0 * M_PI * rho / 3.0;

    for(int bin = 0; bin < MAXBIN; bin++){
        // lower and upper radius 
        double rlower = bin * delr;
        double rupper = rlower + delr;
       // expected number of particles for  ideal gas
        double nideal = constant * (pow(rupper,3) - pow(rlower,3));
        
        // function definition
        // actual pair counts/ideal gas reference
        gr[bin] = hist[bin] / (n_particles * nideal);
    }
}

int main(int argc, char* argv[]){

    assert(packing_fraction > 0.0 && packing_fraction < 1.0);
    assert(diameter > 0.0);
    assert(delta > 0.0);

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    read_data();

    if(n_particles == 0){
        printf("Error: Number of particles, n_particles = 0.\n");
        return 0;
    }

    set_packing_fraction();

    double rmax = box[0] / 2.0;
    delr = rmax / MAXBIN;

    for(int i=0;i<MAXBIN;i++){
        hist[i] = 0.0;
    }

    dsfmt_seed(time(NULL));

    int accepted = 0;
    int step, n;
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            accepted += move_particle();
        }
        gofr();
        if(step % output_steps == 0){
            printf("Step %d. Move acceptance: %lf.\n", step, (double)accepted / (n_particles * output_steps));
            accepted = 0;
            write_data(step);
        }
    }
normalise();

//prints g(r) and r
FILE *fp = fopen("gr.dat","w");

for(int bin=0; bin<MAXBIN; bin++){

    double r = (bin + 0.5) * delr;

    fprintf(fp,"%f %f\n", r, gr[bin]);
}

fclose(fp);
    return 0;
}

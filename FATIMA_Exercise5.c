//#-- ========================================================================================
//#--                Exercise 5 - Bea Fatima, Student ID: 1725181
//#-- ========================================================================================
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 100;

const double packing_fraction = 0.6; 

const double diameter = 1.0; //diameter
const double delta  = 0.05; //try with 0.05 for better acceptance
/* Volume change -deltaV, delta V */
const double deltaV = 2.0; //change in volume
/* Reduced pressure \beta P */ 
const double betaP = 50.0; 
const char* init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0; //particle counter 
double radius; //radius of sphere
double particle_volume; //volume of sphere
double r[N][NDIM]; //coordinates of particles
double box[NDIM]; //box dimensions

/* Functions */
int change_volume(void){
    /*--------- Your code goes here -----------*/
    //#1 (following style of pseudocode on pg 16 of lecture notes)
    double vold = box[0] * box[1] * box[2]; //vol before change

    //old coordinates (used if change is rejected)
    double r_old[N][NDIM];
    for (int i = 0; i < n_particles; i++)
        for (int d = 0; d < NDIM; d++)
            r_old[i][d] = r[i][d];
    
    //create volume change between -deltaV and deltaV
    double volchange = (2.0 * dsfmt_genrand() - 1.0) * deltaV;
    double vnew = vold + volchange; //new vol after change
    //reject negative volume
    if (vnew <= 0.0) {
        return 0;
    }

    double scale = cbrt(vnew / vold); //scale factor, assumption 3d

    //new box size 
    double boxnew[NDIM];
    for (int d=0; d<NDIM; d++){
        boxnew[d] = box[d] * scale;
    }
    //looping over all particles/dimensions
    //scaling coordinates 
    for (int i=0; i<n_particles; i++){
        for (int d=0; d<NDIM; d++){
            r[i][d] *= scale;
        }
    }

for (int i = 0; i < n_particles; i++){
    for (int j = i + 1; j < n_particles; j++){
        double dist2 = 0.0;

        for (int d = 0; d < NDIM; d++){
            double dx = r[i][d] - r[j][d];

            // nearest image convention using the new box size
            if (dx >  0.5 * boxnew[d]) dx -= boxnew[d];
            if (dx < -0.5 * boxnew[d]) dx += boxnew[d];

            dist2 += dx * dx;
        }

        if (dist2 < diameter * diameter){
            // if overlap, reject and restore old coordinates
            for (int ii = 0; ii < n_particles; ii++)
                for (int d = 0; d < NDIM; d++)
                    r[ii][d] = r_old[ii][d];

            return 0;
        }
    }
}
//metropolis acceptance rule
//in this case, energy change is 0 (if there are overlaps, it is nonzero but this is accounted for)
double accrule = -betaP * (vnew - vold) + n_particles * log(vnew / vold);

    if (dsfmt_genrand() < exp(accrule)){
        // accept volume move
        for (int d = 0; d < NDIM; d++)
            box[d] = boxnew[d];
        return 1;
    }

    // reject move and restore old coordinates
    for (int i = 0; i < n_particles; i++)
        for (int d = 0; d < NDIM; d++)
            r[i][d] = r_old[i][d];

    return 0;
}

//code from the previous exercise
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
/*--------- Your code goes here -----------*/
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
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", r[n][d]);
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

int main(int argc, char* argv[]){

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

    dsfmt_seed(time(NULL));
            
    printf("#Step \t Volume \t Move-acceptance\t Volume-acceptance \n");

    int move_accepted = 0;
    int vol_accepted = 0;
    int step, n;
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            move_accepted += move_particle();
        }
        vol_accepted += change_volume();

        if(step % output_steps == 0){
            printf("%d \t %lf \t %lf \t %lf \n", 
                    step, box[0] * box[1] * box[2], 
                    (double)move_accepted / (n_particles * output_steps), 
                    (double)vol_accepted /  output_steps);
            move_accepted = 0;
            vol_accepted = 0;
            write_data(step);
        }
    }

    return 0;
}


#include<math.h>
#include<stdio.h>
#include<stdlib.h>




int number_particles = 216; //total number of particle
double rho = 1.2;
//const double Volume=number_particles/rho;
double l_x = 14.60; //(sqrt(Volume));
double l_y = 14.60;
double dt = 0.005; //time step
double m = 1.e0; //mass of particles
double r_c = 2.e0; //cutoff radius
//Uabp is for Potential Energy , Ke is for kinetic Energy and Fijx,Fijy variables to store Forces
double Ke, Uabp;
double position[216][2];
double velocity[216][2];
double f_now[216][2];
double K_bt = 1.e0;
double T;

void hello_world_from_utiliy(){
    printf("Helllo world from utility.h\n");
}

/***** intialization function will intialize the position and velocity *******/
void intialization() {
    int i = 0;
    double b;
    for (int k = 0; k < 16; k++) {
        for (int j = 0; j < 16; j++) {
            if (i < number_particles) {
                position[i][0] = k / (sqrt(rho));
                position[i][1] = j / (sqrt(rho));
                i++;
            }
        }
    }

    for (int i = 0; i < number_particles; i++) //loop for intializing the value of velocities
    {

        for (int j = 0; j < 2; j++) //loop for dimensionality
        {
            b = drand48();
            velocity[i][j] = (double)(2.e0 * b - 1.e0);
        }
    }
}
/*******Rescaling The Velocity*******/
void rescale_velocity() {
    double vprimex, vprimey, v2, w;
    v2 = 0.e0;
    for (int i = 0; i < number_particles; i++) {

        vprimex = vprimex + velocity[i][0];
        vprimey = vprimey + velocity[i][1];

    }
    vprimex = vprimex / number_particles;
    vprimey = vprimey / number_particles;

    for (int i = 0; i < number_particles; i++) {
        velocity[i][0] = velocity[i][0] - vprimex;
        velocity[i][1] = velocity[i][1] - vprimey;
    }

    for (int i = 0; i < number_particles; i++) {
        v2 = v2 + velocity[i][0] * velocity[i][1];
    }

    v2 = v2/number_particles;
    w = sqrt((2.e0*T)/v2);

    for (int i = 0; i < number_particles; i++) {
        velocity[i][0] = velocity[i][0]*w;
        velocity[i][1] = velocity[i][1]*w;
    }
}

/******Force Calculation*******/
void update_force() {
    double fijx, fijy;
    double f_c = 4.0 * (12.0 / pow(r_c, 13) - 6.0 / pow(r_c, 7));
    double v_c = 4.0 * (1.0 / pow(r_c, 12) - 1.0 / pow(r_c, 6)) + f_c * r_c;
    for (int i = 0; i < number_particles; i++) //loop for number of particels
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        f_now[i][j] = 0; //fn is intialized zero and will be calculated now

    }
    Uabp = 0.0; //potential energy is set zero every time step
    for (int i = 0; i < number_particles - 1; i++) {
        for (int j = i + 1; j < number_particles; j++) {
            double dx = position[i][0] - position[j][0]; //original difference in position in x direction
            double dy = position[i][1] - position[j][1]; //original difference  in position in y direction
            if (fabs(dx) > (l_x / 2.0)) {
                dx = dx - l_x * (dx / fabs(dx));
            }

            if (fabs(dy) > (l_x / 2.0)) {
                dy = dy - l_x * (dy / fabs(dy));
            }

            double rij2 = dx * dx + dy * dy; // distance between i,j pair

            if ((rij2) < (4)) // Minimum cutoff potential condition
            {

                double rij2inv = (1.0 / rij2); //1/r^2
                double rij6 = rij2inv * rij2inv * rij2inv; // (1/r)^6
                fijx = 24.0 * rij6 * (2.0 * rij6 - 1.0) * dx * rij2inv - f_c * dx * sqrt(rij2inv); // Force according to LG in x direction
                fijy = 24.0 * rij6 * (2.0 * rij6 - 1.0) * dy * rij2inv - f_c * dy * sqrt(rij2inv); //Force according to LG in y direction
                f_now[i][0] += fijx; //Force on particle i due to j
                f_now[i][1] += fijy;
                f_now[j][0] += -fijx; //Force on particles j due to i (Newton's 3 law)
                f_now[j][1] += -fijy;
                Uabp += 4. * rij6 * (rij6 - 1.) + f_c * sqrt(rij2) - v_c; //Potential Energy
            }
        }

    }
};

/****Position Update*****/
void update_position() {

    for (int i = 0; i < number_particles; i++) //loop for number of particles
    {
        for (int j = 0; j < 2; ++j) //loop for dimensionality
        {
            position[i][j] += dt * velocity[i][j] + ((dt * dt) / (2. * m)) * f_now[i][j]; //position is upadted using velocity verlet algortihum
            if (position[i][j] > l_x) //PBC
            {
                position[i][j] -= l_x;
            }

            if (position[i][j] < 0) {
                position[i][j] += l_x;
            }
        }
    }
};
/****Velocity Update*****/
void update_velocity() {
    Ke = 0.0; //Intially KE=0 for each time step
    for (int i = 0; i < number_particles; i++) //loop for number of particles
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        velocity[i][j] += (dt / (2.e0)) * (f_now[i][j]); //velocities are upadted using velocity verlet algorithum
        Ke += 0.5 * velocity[i][j] * velocity[i][j]; // Kinetice Energy is upadted
    }
};



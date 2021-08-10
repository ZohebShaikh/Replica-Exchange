#include<math.h>
#include<stdio.h>
#define CURRENT_RAND_MAX 10
int number_particles = 250; //total number of particle
double rho = 1.2;
//const double Volume=number_particles/rho;
double l_x = 14.60; //(sqrt(Volume));
double l_y = 14.60;
double dt = 0.005; //time step
double m = 1.e0; //mass of particles
double r_c = 2.e0; //cutoff radius
//Uabp is for Potential Energy , Ke is for kinetic Energy and Fijx,Fijy variables to store Forces
double Ke1, Uabp1, Ke2, Uabp2, Ke3, Uabp3;
double position1[250][2], position2[250][2], position3[250][2];
double velocity1[250][2], velocity2[250][2], velocity3[250][2];
double f_now1[250][2], f_now2[250][2], f_now3[250][2];
double K_bt = 1.e0;
//double wallT=1.e0;

void hello_world_from_utiliy(){
    printf("Helllo world from utility.h\n");
}

/***** intialization function will intialize the position and velocity *******/
void intialization1() {
    int i = 0;
    for (int k = 0; k < 16; k++) {
        for (int j = 0; j < 16; j++) {
            if (i < number_particles) {
                position1[i][0] = k / (sqrt(rho));
                position1[i][1] = j / (sqrt(rho));
                i++;
            }
        }
    }

    for (int i = 0; i < number_particles; i++) //loop for intializing the value of velocities
    {

        for (int j = 0; j < 2; j++) //loop for dimensionality
        {
            velocity1[i][j] = (double)(2 * (rand() / CURRENT_RAND_MAX) - 1);
        }
    }
}
/*******Rescaling The Velocity*******/
void rescale_velocity1() {
    double vprimex, vprimey;
    for (int i = 0; i < number_particles; i++) {

        vprimex = vprimex + velocity1[i][0];
        vprimey = vprimey + velocity1[i][1];

    }
    vprimex = vprimex / number_particles;
    vprimey = vprimey / number_particles;

    for (int i = 0; i < number_particles; i++) {
        velocity1[i][0] = velocity1[i][0] - vprimex;
        velocity1[i][1] = velocity1[i][1] - vprimey;
    }
}

/******Force Calculation*******/
void update_force1() {
    double fijx, fijy;
    double f_c = 4.0 * (12.0 / pow(r_c, 13) - 6.0 / pow(r_c, 7));
    double v_c = 4.0 * (1.0 / pow(r_c, 12) - 1.0 / pow(r_c, 6)) + f_c * r_c;
    for (int i = 0; i < number_particles; i++) //loop for number of particels
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        f_now1[i][j] = 0; //fn is intialized zero and will be calculated now

    }
    Uabp1 = 0.0; //potential energy is set zero every time step
    for (int i = 0; i < number_particles - 1; i++) {
        for (int j = i + 1; j < number_particles; j++) {
            double dx = position1[i][0] - position1[j][0]; //original difference in position in x direction
            double dy = position1[i][1] - position1[j][1]; //original difference  in position in y direction
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
                f_now1[i][0] += fijx; //Force on particle i due to j
                f_now1[i][1] += fijy;
                f_now1[j][0] += -fijx; //Force on particles j due to i (Newton's 3 law)
                f_now1[j][1] += -fijy;
                Uabp1 += 4. * rij6 * (rij6 - 1.) + f_c * sqrt(rij2) - v_c; //Potential Energy
            }
        }

    }
};

/****Position Update*****/
void update_position1() {

    for (int i = 0; i < number_particles; i++) //loop for number of particles
    {
        for (int j = 0; j < 2; ++j) //loop for dimensionality
        {
            position1[i][j] += dt * velocity1[i][j] + ((dt * dt) / (2. * m)) * f_now1[i][j]; //position is upadted using velocity verlet algortihum
            if (position1[i][j] > l_x) //WALL
            {
                position1[i][j] -= l_x;
            }

            if (position1[i][j] < 0) {
                position1[i][j] += l_x;
            }
        }
    }
};
/****Velocity Update*****/
void update_velocity1() {
    Ke1 = 0.0; //Intially KE=0 for each time step
    for (int i = 0; i < number_particles; i++) //loop for number of particles
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        velocity1[i][j] += (dt / (2)) * (f_now1[i][j]); //velocities are upadted using velocity verlet algorithum
        Ke1 += 0.5 * velocity1[i][j] * velocity1[i][j]; // Kinetice Energy is upadted
    }
};

void intialization2() {
    int i = 0;
    for (int k = 0; k < 16; k++) {
        for (int j = 0; j < 16; j++) {
            if (i < number_particles) {
                position2[i][0] = k / (sqrt(rho));
                position2[i][1] = j / (sqrt(rho));
                i++;
            }
        }
    }

    for (int i = 0; i < number_particles; i++) //loop for intializing the value of velocities
    {

        for (int j = 0; j < 2; j++) //loop for dimensionality
        {
            velocity2[i][j] = (double)(2 * (rand() / CURRENT_RAND_MAX) - 1);
        }
    }
}
/*******Rescaling The Velocity*******/
void rescale_velocity2() {
    double vprimex, vprimey;
    for (int i = 0; i < number_particles; i++) {

        vprimex = vprimex + velocity2[i][0];
        vprimey = vprimey + velocity2[i][1];

    }
    vprimex = vprimex / number_particles;
    vprimey = vprimey / number_particles;

    for (int i = 0; i < number_particles; i++) {
        velocity2[i][0] = velocity2[i][0] - vprimex;
        velocity2[i][1] = velocity2[i][1] - vprimey;
    }
}

/******Force Calculation*******/
void update_force2() {
    double fijx, fijy;
    double f_c = 4.0 * (12.0 / pow(r_c, 13) - 6.0 / pow(r_c, 7));
    double v_c = 4.0 * (1.0 / pow(r_c, 12) - 1.0 / pow(r_c, 6)) + f_c * r_c;
    for (int i = 0; i < number_particles; i++) //loop for number of particels
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        f_now2[i][j] = 0; //fn is intialized zero and will be calculated now

    }
    Uabp2 = 0.0; //potential energy is set zero every time step
    for (int i = 0; i < number_particles - 1; i++) {
        for (int j = i + 1; j < number_particles; j++) {
            double dx = position2[i][0] - position2[j][0]; //original difference in position in x direction
            double dy = position2[i][1] - position2[j][1]; //original difference  in position in y direction
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
                f_now2[i][0] += fijx; //Force on particle i due to j
                f_now2[i][1] += fijy;
                f_now2[j][0] += -fijx; //Force on particles j due to i (Newton's 3 law)
                f_now2[j][1] += -fijy;
                Uabp2 += 4. * rij6 * (rij6 - 1.) + f_c * sqrt(rij2) - v_c; //Potential Energy
            }
        }

    }
};

/****Position Update*****/
void update_position2() {

    for (int i = 0; i < number_particles; i++) //loop for number of particles
    {
        for (int j = 0; j < 2; ++j) //loop for dimensionality
        {
            position2[i][j] += dt * velocity2[i][j] + ((dt * dt) / (2. * m)) * f_now2[i][j]; //position is upadted using velocity verlet algortihum
            if (position2[i][j] > l_x) //WALL
            {
                position2[i][j] -= l_x;
            }

            if (position2[i][j] < 0) {
                position2[i][j] += l_x;
            }
        }
    }
};
/****Velocity Update*****/
void update_velocity2() {
    Ke2 = 0.0; //Intially KE=0 for each time step
    for (int i = 0; i < number_particles; i++) //loop for number of particles
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        velocity2[i][j] += (dt / (2)) * (f_now2[i][j]); //velocities are upadted using velocity verlet algorithum
        Ke2 += 0.5 * velocity2[i][j] * velocity2[i][j]; // Kinetice Energy is upadted
    }
};

void intialization3() {
    int i = 0;
    for (int k = 0; k < 16; k++) {
        for (int j = 0; j < 16; j++) {
            if (i < number_particles) {
                position3[i][0] = k / (sqrt(rho));
                position3[i][1] = j / (sqrt(rho));
                i++;
            }
        }
    }

    for (int i = 0; i < number_particles; i++) //loop for intializing the value of velocities
    {

        for (int j = 0; j < 2; j++) //loop for dimensionality
        {
            velocity3[i][j] = (double)(2 * (rand() / CURRENT_RAND_MAX) - 1);
        }
    }
}
/*******Rescaling The Velocity*******/
void rescale_velocity3() {
    double vprimex, vprimey;
    for (int i = 0; i < number_particles; i++) {

        vprimex = vprimex + velocity3[i][0];
        vprimey = vprimey + velocity3[i][1];

    }
    vprimex = vprimex / number_particles;
    vprimey = vprimey / number_particles;

    for (int i = 0; i < number_particles; i++) {
        velocity3[i][0] = velocity3[i][0] - vprimex;
        velocity3[i][1] = velocity3[i][1] - vprimey;
    }
}

/******Force Calculation*******/
void update_force3() {
    double fijx, fijy;
    double f_c = 4.0 * (12.0 / pow(r_c, 13) - 6.0 / pow(r_c, 7));
    double v_c = 4.0 * (1.0 / pow(r_c, 12) - 1.0 / pow(r_c, 6)) + f_c * r_c;
    for (int i = 0; i < number_particles; i++) //loop for number of particels
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        f_now3[i][j] = 0; //fn is intialized zero and will be calculated now

    }
    Uabp3 = 0.0; //potential energy is set zero every time step
    for (int i = 0; i < number_particles - 1; i++) {
        for (int j = i + 1; j < number_particles; j++) {
            double dx = position3[i][0] - position3[j][0]; //original difference in position in x direction
            double dy = position3[i][1] - position3[j][1]; //original difference  in position in y direction
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
                f_now3[i][0] += fijx; //Force on particle i due to j
                f_now3[i][1] += fijy;
                f_now3[j][0] += -fijx; //Force on particles j due to i (Newton's 3 law)
                f_now3[j][1] += -fijy;
                Uabp3 += 4. * rij6 * (rij6 - 1.) + f_c * sqrt(rij2) - v_c; //Potential Energy
            }
        }

    }
};

/****Position Update*****/
void update_position3() {

    for (int i = 0; i < number_particles; i++) //loop for number of particles
    {
        for (int j = 0; j < 2; ++j) //loop for dimensionality
        {
            position3[i][j] += dt * velocity3[i][j] + ((dt * dt) / (2. * m)) * f_now3[i][j]; //position is upadted using velocity verlet algortihum
            if (position3[i][j] > l_x) //WALL
            {
                position3[i][j] -= l_x;
            }

            if (position3[i][j] < 0) {
                position3[i][j] += l_x;
            }
        }
    }
};
/****Velocity Update*****/
void update_velocity3() {
    Ke3 = 0.0; //Intially KE=0 for each time step
    for (int i = 0; i < number_particles; i++) //loop for number of particles
        for (int j = 0; j < 2; j++) //loop for dimensionality
    {
        velocity3[i][j] += (dt / (2)) * (f_now3[i][j]); //velocities are upadted using velocity verlet algorithum
        Ke3 += 0.5 * velocity3[i][j] * velocity3[i][j]; // Kinetice Energy is upadted
    }
};
#include <stdio.h>

#include <stdlib.h>

#include "utility.h"

#include <mpi.h>

#define MAX_ITERATIONS 10001
#define MASTER_RANK 0
#define WORKER_1_RANK 1
#define WORKER_2_RANK 2
#define WORKER_3_RANK 3

int main(int argc, char ** argv) {

    // Unique rank is assigned to each process in a communicator
    int rank;

    // Total number of ranks
    int size;

    // Initializes the MPI execution environment
    MPI_Init( & argc, & argv);

    // Get this process' rank (process within a communicator)
    // MPI_COMM_WORLD is the default communicator
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);

    // Get the total number ranks in this communicator
    MPI_Comm_size(MPI_COMM_WORLD, & size);

    // Main Logic start
    double * temperature_1 = (double * ) malloc(sizeof(double));
    double * temperature_2 = (double * ) malloc(sizeof(double));
    double * temperature_3 = (double * ) malloc(sizeof(double));
    int temp = 0;
    double t;
    double de; //difference in energy
    double dBeta; //difference in inverses of temperature
    // Synchronize so we can remove interleaved output
    if (rank == 0) {
        double position1[216][2], position2[216][2], position3[216][2];
        double velocity1[216][2], velocity2[216][2], velocity3[216][2];
        double E1, E2, E3;
        int count;
        // Master Logic
        // Always print from rank 0
        while (1) {
            //if (count == MAX_ITERATIONS){
            // Run till count is less than max iterations
            MPI_Recv( & count, 1, MPI_INT, WORKER_1_RANK, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Count:  %d\n", count); // try to playaround with this break don't know why count is not working properly
            if (count == 10000) break;
            // Takes buffer, size, type, source, tag, communicator, and status

            printf("Waiting for 1, 2 and 3\n");
            // //         int MPI_Recv(void *buf, int count, MPI_Datatype datatype,
            // // int source, int tag, MPI_Comm comm, MPI_Status *status)
            MPI_Recv(temperature_1, 1, MPI_DOUBLE, WORKER_1_RANK, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //MPI_Rec(the energy ("Uabp") from rank 1 and put in "E1")
            MPI_Recv( & E1, 1, MPI_DOUBLE, WORKER_1_RANK, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("The energy 1 is %lf\n", E1);

            //MPI_Rec(the position array from rank 1 and put in "position1")
            MPI_Recv( & position1, 216 * 2, MPI_DOUBLE, WORKER_1_RANK, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //MPI_Rec(the velocity array from rank 1 and put in "velocity1")
            MPI_Recv( & velocity1, 216 * 2, MPI_DOUBLE, WORKER_1_RANK, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            printf("Waiting for 2 and 3\n");
            MPI_Recv(temperature_2, 1, MPI_DOUBLE, WORKER_2_RANK, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //MPI_Rec(the energy ("Uabp") from rank 2 and put in "E2")
            MPI_Recv( & E2, 1, MPI_DOUBLE, WORKER_2_RANK, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //MPI_Rec(the position array from rank 2 and put in "position2")
            MPI_Recv( & position2, 216 * 2, MPI_DOUBLE, WORKER_2_RANK, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //MPI_Rec(the velocity array from rank 2 and put in "velocity2")
            MPI_Recv( & velocity2, 216 * 2, MPI_DOUBLE, WORKER_2_RANK, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            printf("Waiting for 3\n");
            MPI_Recv(temperature_3, 1, MPI_DOUBLE, WORKER_3_RANK, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            printf("Waiting complete\n");

            //MPI_Rec(the energy ("Uabp") from rank 3 and put in "E3")
            MPI_Recv( & E3, 1, MPI_DOUBLE, WORKER_3_RANK, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //MPI_Rec(the position array from rank 3 and put in "position3")
            MPI_Recv( & position3, 216 * 2, MPI_DOUBLE, WORKER_3_RANK, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            //MPI_Rec(the velocity array from rank 3 and put in "velocity3")
            MPI_Recv( & velocity3, 216 * 2, MPI_DOUBLE, WORKER_3_RANK, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Print our received message
            printf("%f\n", * temperature_1);
            printf("%f\n", * temperature_2);
            printf("%f\n", * temperature_3);
            //printf("Count:  %d\n",count);
            de = 0.e0;
            dBeta = 0.e0;

            if ((count % 2000) != 0) //  1000,3000
            {
                de = E1 - E2;
                dBeta = 1.e0 / ( * temperature_1) - 1.e0 / ( * temperature_2);
                if (drand48() < exp(-de * dBeta)) {

                    for (int i = 0; i < number_particles; i++) //loop for exchange postions
                    {
                        for (int j = 0; j < 2; j++) //loop for dimensionality
                        {
                            velocity1[i][j] = sqrt( * temperature_2 / ( * temperature_1)) * velocity1[i][j];
                            velocity2[i][j] = sqrt(( * temperature_1) / ( * temperature_2)) * velocity2[i][j];
                        }
                    }

                    //  MPI_SEND(&position2,,,,,to WORKER_1_RANK)
                    MPI_Send( & position2, 216 * 2, MPI_DOUBLE, WORKER_1_RANK, 13, MPI_COMM_WORLD);
                    //    MPI_SEND(&velocity2,,,,,to WORKER_1_RANK);
                    MPI_Send( & velocity2, 216 * 2, MPI_DOUBLE, WORKER_1_RANK, 14, MPI_COMM_WORLD);
                    //  MPI_SEND(&position1,,,,,to WORKER_2_RANK)
                    MPI_Send( & position1, 216 * 2, MPI_DOUBLE, WORKER_2_RANK, 15, MPI_COMM_WORLD);
                    //    MPI_SEND(&velocity1,,,,,to WORKER_2_RANK);
                    MPI_Send( & velocity1, 216 * 2, MPI_DOUBLE, WORKER_2_RANK, 16, MPI_COMM_WORLD);
                    //  MPI_SEND(&position3,,,,,to WORKER_3_RANK) (No swap for Replica 3)
                    MPI_Send( & position3, 216 * 2, MPI_DOUBLE, WORKER_3_RANK, 17, MPI_COMM_WORLD);
                    //    MPI_SEND(&velocity3,,,,,to WORKER_3_RANK);
                    MPI_Send( & velocity3, 216 * 2, MPI_DOUBLE, WORKER_3_RANK, 18, MPI_COMM_WORLD);
                }
            } else { // 2000 , 4000, 6000
                de = E2 - E3;
                dBeta = 1.e0 / ( * temperature_2) - 1.e0 / ( * temperature_3);
                if (drand48() < exp(-de * dBeta)) {

                    for (int i = 0; i < number_particles; i++) //loop for exchange postions
                    {
                        for (int j = 0; j < 2; j++) //loop for dimensionality
                        {
                            velocity2[i][j] = sqrt(( * temperature_3) / ( * temperature_2)) * velocity2[i][j];
                            velocity3[i][j] = sqrt(( * temperature_2) / ( * temperature_3)) * velocity3[i][j];
                        }
                    }

                    //  MPI_SEND(&position1,,,,,to WORKER_1_RANK) (No Swap here for Replica 1)
                    MPI_Send( & position1, 216 * 2, MPI_DOUBLE, WORKER_1_RANK, 13, MPI_COMM_WORLD);
                    //    MPI_SEND(&velocity1,,,,,to WORKER_1_RANK);
                    MPI_Send( & velocity1, 216 * 2, MPI_DOUBLE, WORKER_1_RANK, 14, MPI_COMM_WORLD);
                    // MPI_SEND(&position3,,,,,to WORKER_2_RANK);
                    MPI_Send( & position3, 216 * 2, MPI_DOUBLE, WORKER_2_RANK, 15, MPI_COMM_WORLD);
                    // MPI_SEND(&velocity3,,,,,to WORKER_2_RANK);
                    MPI_Send( & velocity3, 216 * 2, MPI_DOUBLE, WORKER_2_RANK, 16, MPI_COMM_WORLD);
                    // MPI_SEND(&position2,,,,,to WORKER_3_RANK);
                    MPI_Send( & position2, 216 * 2, MPI_DOUBLE, WORKER_3_RANK, 17, MPI_COMM_WORLD);
                    // MPI_SEND(&velocity2,,,,,to WORKER_3_RANK);
                    MPI_Send( & velocity2, 216 * 2, MPI_DOUBLE, WORKER_3_RANK, 18, MPI_COMM_WORLD);
                }
            }
            //MPI_Barrier(MPI_COMM_WORLD);
            //break; // try to playaround with this break don't know why count is not working properly
        }
    }
    // Replica 1
    else if (rank == 1) {
        double pos[216][2];
        double vel[216][2];
        FILE * fp11 = fopen("KE1.txt", "w");
        FILE * fp21 = fopen("PE1.txt", "w");
        FILE * fp31 = fopen("TE1.txt", "w");
        int WORKER_count;
        WORKER_count = 0;
        T = 0.5e0;
        intialization();
        rescale_velocity();
        update_force();
        for (int i = 0; i < MAX_ITERATIONS; i++) {
            t=i*dt;
            update_position();
            update_velocity();
            update_force();
            update_velocity();
            if ((i % 1000) == 0) {
                //sending counter
                MPI_Send( & WORKER_count, 1, MPI_INT, MASTER_RANK, 0, MPI_COMM_WORLD);
                //sending temperature
                MPI_Send( & T, 1, MPI_DOUBLE, MASTER_RANK, 1, MPI_COMM_WORLD);
                //MPI_Send(the value of energy);
                MPI_Send( & Uabp, 1, MPI_DOUBLE, MASTER_RANK, 2, MPI_COMM_WORLD);
                //MPI_Send(the position array);
                MPI_Send( & position, 216 * 2, MPI_DOUBLE, MASTER_RANK, 3, MPI_COMM_WORLD);
                //MPI_Send(the velocity array);
                MPI_Send( & velocity, 216 * 2, MPI_DOUBLE, MASTER_RANK, 4, MPI_COMM_WORLD);
                //MPI_Rec(the new positions);
                MPI_Recv( & pos, 216 * 2, MPI_DOUBLE, MASTER_RANK, 13, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //MPI_Rec(the new velocities);
                MPI_Recv( & vel, 216 * 2, MPI_DOUBLE, MASTER_RANK, 14, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            WORKER_count++;
            if (i % 100 == 0) {
                fprintf(fp11, "%lf\t%lf\n", t, Ke / number_particles);
                fprintf(fp21, "%lf\t%lf\n", t, Uabp / number_particles);
                fprintf(fp31, "%lf\t%lf\n", t, (Ke + Uabp) / number_particles);
            }
        }
        fclose(fp11);
        fclose(fp21);
        fclose(fp31);
    }

    // Replica 2
    else if (rank == 2) {
        double pos[216][2];
        double vel[216][2];
        FILE * fp12 = fopen("KE2.txt", "w");
        FILE * fp22 = fopen("PE2.txt", "w");
        FILE * fp32 = fopen("TE2.txt", "w");

        T = 0.6e0;
        intialization();
        rescale_velocity();
        update_force();
        for (int i = 0; i < MAX_ITERATIONS; i++) {
            t = i * dt;
            update_position();
            update_velocity();
            update_force();
            update_velocity();
            if ((i % 1000) == 0) {
                //sending temperature
                MPI_Send( & T, 1, MPI_DOUBLE, MASTER_RANK, 5, MPI_COMM_WORLD);
                //MPI_Send(the value of energy);
                MPI_Send( & Uabp, 1, MPI_DOUBLE, MASTER_RANK, 6, MPI_COMM_WORLD);
                //MPI_Send(the position array);
                MPI_Send( & position, 216 * 2, MPI_DOUBLE, MASTER_RANK, 7, MPI_COMM_WORLD);
                //MPI_Send(the velocity array);
                MPI_Send( & velocity, 216 * 2, MPI_DOUBLE, MASTER_RANK, 8, MPI_COMM_WORLD);
                //MPI_Rec(the new positions);
                MPI_Recv( & pos, 216 * 2, MPI_DOUBLE, MASTER_RANK, 15, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //MPI_Rec(the new velocities);
                MPI_Recv( & vel, 216 * 2, MPI_DOUBLE, MASTER_RANK, 16, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            if (i % 100 == 0) {
                fprintf(fp12, "%lf\t%lf\n", t, Ke / number_particles);
                fprintf(fp22, "%lf\t%lf\n", t, Uabp / number_particles);
                fprintf(fp32, "%lf\t%lf\n", t, (Ke + Uabp) / number_particles);
            }
        }
        fclose(fp12);
        fclose(fp22);
        fclose(fp32);
    }

    // Replica 3
    else if (rank == 3) {
        double pos[216][2];
        double vel[216][2];
        FILE * fp13 = fopen("KE3.txt", "w");
        FILE * fp23 = fopen("PE3.txt", "w");
        FILE * fp33 = fopen("TE3.txt", "w");

        T = 0.7e0;
        intialization();
        rescale_velocity();
        update_force();
        for (int i = 0; i < MAX_ITERATIONS; i++) {
            t = i * dt;
            update_position();
            update_velocity();
            update_force();
            update_velocity();
            if ((i % 1000) == 0) {
                //sending counter
                // MPI_Send( WORKER_count , 1, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD);
                //sending temperature
                MPI_Send( & T, 1, MPI_DOUBLE, MASTER_RANK, 9, MPI_COMM_WORLD);
                //MPI_Send(the value of energy);
                MPI_Send( & Uabp, 1, MPI_DOUBLE, MASTER_RANK, 10, MPI_COMM_WORLD);
                //MPI_Send(the position array);
                MPI_Send( & position, 216 * 2, MPI_DOUBLE, MASTER_RANK, 11, MPI_COMM_WORLD);
                //MPI_Send(the velocity array);
                MPI_Send( & velocity, 216 * 2, MPI_DOUBLE, MASTER_RANK, 12, MPI_COMM_WORLD);
                //MPI_Rec(the new positions);
                MPI_Recv( & pos, 216 * 2, MPI_DOUBLE, MASTER_RANK, 17, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //MPI_Rec(the new velocities);
                MPI_Recv( & vel, 216 * 2, MPI_DOUBLE, MASTER_RANK, 18, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (i % 100 == 0) {
                fprintf(fp13, "%lf\t%lf\n", t, Ke / number_particles);
                fprintf(fp23, "%lf\t%lf\n", t, Uabp / number_particles);
                fprintf(fp33, "%lf\t%lf\n", t, (Ke + Uabp) / number_particles);
            }
        }
        fclose(fp13);
        fclose(fp23);
        fclose(fp33);
    }
    // Terminate MPI execution environment
    MPI_Finalize();
}
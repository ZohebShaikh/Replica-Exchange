#include <stdio.h>
#include <stdlib.h>
#include "utility.h"
#include <mpi.h>

#define MAX_ITERATIONS 10000
#define MASTER_RANK 0
#define SLAVE_1_RANK 1
#define SLAVE_2_RANK 2
#define SLAVE_3_RANK 3

int main(int argc, char ** argv) {
    
    // Checking utility working
    hello_world_from_utiliy();

    // Unique rank is assigned to each process in a communicator
    int rank;

    // Total number of ranks
    int size;

    // The machine we are on
    char name[80];

    // Length of the machine name
    int length;

    // Initializes the MPI execution environment
    MPI_Init( & argc, & argv);

    // Get this process' rank (process within a communicator)
    // MPI_COMM_WORLD is the default communicator
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);

    // Get the total number ranks in this communicator
    MPI_Comm_size(MPI_COMM_WORLD, & size);

    // Gets the name of the processor
    // Implementation specific (may be gethostname, uname, or sysinfo)
    MPI_Get_processor_name(name, & length);

    // Main Logic start
    float *temperature_1 = (float *) malloc ( sizeof ( float) ) ;
    float *temperature_2 = (float *) malloc ( sizeof ( float) ) ;
    float *temperature_3 = (float *) malloc ( sizeof ( float) ) ;
    int temp = 0;
    int *count = (int * ) malloc (sizeof(int)); 
    *count = 0;
    // Synchronize so we can remove interleaved output
    if (rank == 0) {
        // Master Logic
        // Always print from rank 0
        while (1){
            // Run till count is less than max iterations
            MPI_Recv(count, 1, MPI_INT, SLAVE_1_RANK, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //if (*count == MAX_ITERATIONS) break;
            printf("Count:  %d\n",*count);  // try to playaround with this break don't know why count is not working properly
            // Takes buffer, size, type, source, tag, communicator, and status
            printf("Waiting for 1, 2 and 3\n");
            MPI_Recv(temperature_1, 1, MPI_FLOAT, SLAVE_1_RANK, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Waiting for 2 and 3\n");
            MPI_Recv(temperature_2, 1, MPI_FLOAT, SLAVE_2_RANK, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Waiting for 3\n");
            MPI_Recv(temperature_3, 1, MPI_FLOAT, SLAVE_3_RANK, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Waiting complete\n");
            // Print our received message
            printf("%f\n", *temperature_1);
            printf("%f\n", *temperature_2);
            printf("%f\n", *temperature_3);

            //temp  = count / 1000;

            //if ((count%2) == 0  ) // even 2000,40000{
                // swap(T2,T3);
                // MPI Send T2
            //}   
            //else{ // 1000 , 3000,5000
                // swap(T1,T2);
                // MPI send T1
                // 
            //}
            break; // try to playaround with this break don't know why count is not working properly
        }
    } else {
        // If not rank zero, work as slave
        float temperature = 12.0 * rank;
        int *slave_count = (int * ) malloc (sizeof(int));
        *slave_count = 0;
        for ( int i = 0 ; i < MAX_ITERATIONS; i++){
            *slave_count ++ ;
            if ( (i % 1000 )== 0){
                    MPI_Send( slave_count , 1, MPI_INT, MASTER_RANK, rank, MPI_COMM_WORLD);
                    MPI_Send( &temperature , 1, MPI_FLOAT, MASTER_RANK, rank, MPI_COMM_WORLD);
            }
        }
    }
    // Terminate MPI execution environment
    MPI_Finalize();
}
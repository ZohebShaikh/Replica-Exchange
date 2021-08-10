#include <stdio.h>
#include<string.h>
#include <mpi.h>

int temp_1 = 0;
int temp_2 = 0;
int temp_3 = 0;

int main(int argc, char ** argv) {

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

    // Pack these values together into a string
    int buffer_len = 150;
    char buffer[buffer_len];
    sprintf(buffer, "Hello, MPI! Rank: %d Total: %d Machine: %s", rank, size,
        name);
    char temp_1[buffer_len];
    char temp_2[buffer_len];
    char temp_3[buffer_len];
    memset(temp_1, '\0', sizeof(temp_1));
    memset(temp_2, '\0', sizeof(temp_2));
    memset(temp_3, '\0', sizeof(temp_3));

    // Synchronize so we can remove interleaved output
    if (rank == 0) {
        // Always print from rank 0
        puts(buffer);
        printf("Size : %d\n",size);
            // Takes buffer, size, type, source, tag, communicator, and status
            MPI_Recv(buffer, buffer_len, MPI_CHAR, 1, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            strcpy(temp_1,buffer);
            printf("Waiting for 2 and 3\n");
            MPI_Recv(buffer, buffer_len, MPI_CHAR, 2, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            strcpy(temp_2,buffer);
            printf("Waiting for 3\n");
            MPI_Recv(buffer, buffer_len, MPI_CHAR, 3, MPI_ANY_TAG,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("Waiting complete\n");
            strcpy(temp_3,buffer);
            // Print our received message
            printf("%s\n", temp_1);
            printf("%s\n", temp_2);
            printf("%s\n", temp_3);
    } else {
        // If not rank zero, send your message to be printed
        int temp;
        temp =rank;
        size = temp*10;
        for ( int i = 0 ; i < 10000000; i++){
            ;//do some work
        }
        sprintf(buffer, "Hello, MPI! Rank: %d Temp: %d Machine: %s", rank, size,
        name);
        
        //MPI_Send(buffer,buffer_len,datatype,destination,tag,comms)
        MPI_Send(buffer, buffer_len, MPI_CHAR, 0, rank, MPI_COMM_WORLD);
    }

    // Terminate MPI execution environment
    MPI_Finalize();
}
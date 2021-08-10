    MPI_Init(&argc, &argv);

    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0){
        int x;
        MPI_Send(&x, MPI_INT,1,msgtag,
        MPI_COMM_WORLD);
        master();

    }
    else{        
        int x;
        //msg tag can be MPI_ANY_TAG
        //synchronosz send/recive
        //MPI Send
        // MPI_Send(buf,count,datatype,dest,tag,comm=MPI_COMM_WORLD);
        MPI_Send(&x, MPI_INT,0,msgtag,
        MPI_COMM_WORLD);
        slave();
    }
    MPI_Finalize();
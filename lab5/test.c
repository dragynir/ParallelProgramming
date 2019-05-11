#include <pthread.h>
#include <stdio.h>
#include<stdlib.h>
#include<mpi.h>


int main(int argc, char** argv)
{
    MPI_Init(&argc  , &argv);
    int m_rank, p_count;
    MPI_Comm_size(MPI_COMM_WORLD , &p_count);
    MPI_Comm_rank(MPI_COMM_WORLD , &m_rank);

    //if(m_rank == 0)return 0;


    int count = 1;
    MPI_Status status;
    //printf("sdfsdfdsfdsf\n");
    int code = 9999;
    /*if(m_rank == 1)
        code = MPI_Recv(&count, 1, MPI_INT, 0, 12, MPI_COMM_WORLD, &status);*/

    if(m_rank == 0){
        code = MPI_Send(&count, 1, MPI_INT, 1, 12, MPI_COMM_WORLD);
        printf("AA%d\n" , code);

    }
    printf("%d\n" , code);

    


    //int code =  MPI_Send(&count, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);


    //if(code != MPI_SUCCESS)printf("AAAAAa\n");
    //printf("%d\n", code);
    MPI_Finalize();

    return 0;
}


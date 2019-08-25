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


    if(m_rank ==0 ){
        return 0;
    }

    printf("%d\n" , m_rank);

    
    MPI_Finalize();

    return 0;
}


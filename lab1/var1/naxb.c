#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>


//58 : 58


void dotProduct(double* a , double* b , int size , double* dest){
	double s = 0;
	for(int i = 0; i < size; ++i){
		s+=(a[i] * b[i]);
	}
	*dest = s;
}


void absoluteValue(double *a , int size ,  double* dest){
	double d = 0;
	dotProduct(a , a , size , &d);
	*dest = sqrt(d);
}


void multiMatVec(double* A ,  double* x , double* dest , int A_rows , int A_col){
	for(int r = 0; r < A_rows; ++r){
		double s = 0;
		for(int c = 0; c < A_col; ++c){
			s+=(A[r * A_col + c] * x[c]); 
		}
		dest[r] = s;
	}
}

void subtractVec(double* a , double* b , double* dest , int size){
	for(int i = 0; i < size; ++i){
		dest[i] = (a[i] - b[i]);
	}
}

void addVec(double* a , double* b , double* dest , int size){
	for(int i = 0; i < size; ++i){
		dest[i] = (a[i] + b[i]);
	}
}


void fillVector(double* v , int size , double value){
	for (int i = 0; i < size; ++i)
	{
		v[i] = value;
	}
}


void multiOnScalar(double* v , double scalar , int size){
	for (int i = 0; i < size; ++i)
	{
		v[i]*=scalar;
	}
}

void showVector(double* v , int size){
	for (int i = 0; i < size; ++i)
	{
		printf("%lf\n", v[i]);
	}
}


double checkAccuracy(double* Ax , double* b , int size){
	double a1 = 0 , a2 = 0;
	absoluteValue(Ax , size , &a1);
	absoluteValue(b , size , &a2);
	return a1 / a2;
}

void fillVectorM3(double* v , int size){
	for (int i = 0; i < size; ++i)
	{
		v[i] = 1;
		if(!(i%3)) v[i] = (1 - 2*(i%2))*50; 
		else v[i] = 0;
	}
}

void fillVectorM3Grid(double* grid , double* b , int N_x , int N_y){
	for (int i = 0; i < N_y; ++i)
	{
		for (int j = 0; j < N_x; ++j)
		{
			grid[i * N_x + j] = b[i * N_x + j];
		}
	}
}

void fillMatrixM3(double* m , int r , int c , int l_s_index , int N_x){



	for(int i = 0; i < r; ++i){
		for(int j = 0; j < c; ++j){
			if(j == l_s_index + i){
				m[i * c + j] = -4;
			}else if(j == l_s_index + i + N_x || j == l_s_index + i - N_x){
				m[i * c + j] = 1;
			}else if(j == l_s_index	+ i + 1){
				if(j && (j % N_x == 0)){
					m[i * c + j] = 0;
				}else{
					m[i * c + j] = 1;
				}
			}else if(j == l_s_index	+ i - 1){
				if(j && ((j + 1) % N_x == 0)){
					m[i * c + j] = 0;
				}else{
					m[i * c + j] = 1;
				}
			}
			else{
				m[i * c + j] = 0;
			}
		}
	}
}


void showMatrix(double* m , int r  , int c){
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			printf("%d ", (int)m[i * c + j]);
		}
		printf("\n");
	}
}



int main(int argc , char** argv){
	MPI_Init(&argc  , &argv);
	srand(time(0));
	if(argc != 3){
		printf("N e a\n");
		return 0;
	}




	int N_x = atoi(argv[1]);
	int N_y = atoi(argv[2]);
	if(N_x <= 0 || N_y <= 0){
		 printf("Wrong Args\n");
		 return 0;
	}

	int MATRIX_SIZE = N_x * N_y;

	int m_rank = 0 , p_count = 0;

	MPI_Comm_size(MPI_COMM_WORLD , &p_count);
	MPI_Comm_rank(MPI_COMM_WORLD , &m_rank);

	if(0 == m_rank){
		printf("Process count: %d\n", p_count);
	}


	int task_per_process = MATRIX_SIZE / p_count;
	int overflow_tasks = MATRIX_SIZE % p_count;
	int l_s_index = 0 , l_e_index = 0;
	l_s_index = m_rank * task_per_process;
	l_e_index = l_s_index + task_per_process - 1;


	int* recvcounts = (int*)malloc(sizeof(int) * p_count);
	int* displs = (int*)malloc(sizeof(int) * p_count);
	int i_sum = 0;

	for (int i = 0; i < p_count; ++i)
	{
		double s = i * task_per_process + (i < overflow_tasks ? i : overflow_tasks);
		double e = s + task_per_process - 1 + (i < overflow_tasks ? 1 : 0);
		displs[i] = i_sum;
		recvcounts[i] = e - s + 1;
		i_sum+=recvcounts[i];
		if(m_rank == i){
			l_s_index = s;l_e_index = e;
		}
	}

	int l_matrix_r = (l_e_index - l_s_index + 1);
	int l_matrix_c = MATRIX_SIZE;
	double* l_matrix = (double*)malloc(l_matrix_c * l_matrix_r * sizeof(double));
	double* l_vec_keepAx = (double*)malloc(sizeof(double) * l_matrix_r);
	double* vec_x = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_b = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_Gathered_y = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_Gathered = (double*)malloc(sizeof(double) * MATRIX_SIZE);	
	double  t_factor = 0;


	/////////TempTest
	fillVectorM3(vec_b , MATRIX_SIZE);
	fillMatrixM3(l_matrix , l_matrix_r , l_matrix_c , l_s_index , N_x);
	fillVector(vec_x , MATRIX_SIZE , 0);
	

	double start_t = MPI_Wtime() , end_t = 0;
	int end = 0;
	int count = 0;
	while(1){
		++count;
		multiMatVec(l_matrix ,  vec_x , l_vec_keepAx , l_matrix_r , l_matrix_c);
		//Reduce
		MPI_Allgatherv(l_vec_keepAx , l_matrix_r , MPI_DOUBLE,vec_Gathered ,
		recvcounts , displs ,  MPI_DOUBLE, MPI_COMM_WORLD);

		//calculate y(n)
		subtractVec(vec_Gathered , vec_b , vec_Gathered_y , MATRIX_SIZE);
		double check = 1;
		if(0 == m_rank){
			check = checkAccuracy(vec_Gathered_y , vec_b , MATRIX_SIZE);
			if(check < 0.00001){
				end = 1;
				printf("res: %lf\n", check);
			}
		}
		MPI_Bcast(&end , 1 , MPI_INT, 0 , MPI_COMM_WORLD);
		
		if(end){
			end_t = MPI_Wtime();	
			if(0 == m_rank){
				//showVector(vec_x , MATRIX_SIZE);
				printf("Time taken: %f  , Iterations: %d \n", end_t - start_t , count);
			}
			
			{
			free(l_matrix);free(l_vec_keepAx);
			free(vec_x);free(vec_b);
			free(vec_Gathered);free(vec_Gathered_y);
			}

			MPI_Finalize();
			return 0;
		}

		//calculate A * y(n)
		multiMatVec(l_matrix , vec_Gathered_y , l_vec_keepAx , l_matrix_r , l_matrix_c);
		
		MPI_Allgatherv(l_vec_keepAx , l_matrix_r , MPI_DOUBLE,vec_Gathered ,
		recvcounts , displs ,  MPI_DOUBLE, MPI_COMM_WORLD);

		double dotpr = 0;
		dotProduct(vec_Gathered_y , vec_Gathered , MATRIX_SIZE , &dotpr);
		t_factor = dotpr;
		dotProduct(vec_Gathered , vec_Gathered , MATRIX_SIZE , &dotpr);
		t_factor/=dotpr;

		multiOnScalar(vec_Gathered_y , t_factor , MATRIX_SIZE);

		subtractVec(vec_x , vec_Gathered_y , vec_x , MATRIX_SIZE);
	}

	return 0;

}


#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
//#define MATRIX_SIZE 15000


void multiMatVecAdd(double* A ,  double* x , double* dest , int A_rows , int A_col , int r_diff , int c_diff){
	for(int r = 0; r < A_rows; ++r){
		double s = 0;
		for(int c = 0; c < A_col; ++c){
			s+=(A[r * r_diff + c_diff + c] * x[c]); 
		}
		dest[r]+=s;
	}
}

void subtractVec(double* a , double* b , double* dest , int size){
	for(int i = 0; i < size; ++i){
		dest[i] = (a[i] - b[i]);
	}
}


void multiOnScalar(double* v , double scalar , int size){
	for (int i = 0; i < size; ++i)
	{
		v[i]*=scalar;
	}
}

void fillMatrixM1(double* m , int r , int c , int l_s_index){

	for(int i = 0; i < r; ++i){
		for(int j = 0; j < c; ++j){
			if((j == l_s_index + i)){
				m[i * c + j] = 2;
			}else{
				m[i * c + j] = 1;
			}
		}
	}
}


void fillVector(double* v , int size , double value){
	for (int i = 0; i < size; ++i)
	{
		v[i] = value;
	}
}


void showMatrix(double* A , int r , int c){
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			printf("%lf ",	A[c * i + j]);
		}
		printf("\n");
	}
}

void showVector(double* v , int size){
	for (int i = 0; i < size; ++i)
	{
		printf("%lf\n", v[i]);
	}
}


void dotProduct(double* a , double* b , int size , double* dest , int ra){
	double s = 0;
	for(int i = 0; i < size; ++i){
		s+=(a[i] * b[i]);
	}
	*dest = s;
}


void copyVector(double* a , double* b , int size){
	for (int i = 0; i < size; ++i)
	{
		a[i] = b[i];
	}
}



void fillVectorM2(double* v , int size){
	for (int i = 0; i < size; ++i)
	{
		//v[i] = 0.05;
		v[i] = sin((3.14 * 2 * i) / size);
	}
}





void subToTest(double* a , double* b , double* dest , int size , int b_start){
	for (int i = 0; i < size; ++i)
	{
		dest[i] = (a[i] - b[i + b_start]);
	}
}


void multiMatVecTest(double* A ,  double* x , double* dest , int A_rows , int A_col){
	for(int r = 0; r < A_rows; ++r){
		double s = 0;
		for(int c = 0; c < A_col; ++c){
			s+=(A[r * A_col + c] * x[c]); 
		}
		dest[r] = s;
	}
}



void multi(double* l_matrix , double* l_vec_ , double* l_vec_dest , double* l_vec_recv ,\
 int* recvcounts , int* displs ,   int l_matrix_r , int l_matrix_c ,  int m_rank , int p_count , int task_p_p){
	
	double* l_vec_repr = l_vec_;
	double* l_vec_recv_repr = l_vec_recv;
	int block_i = 0 , exc_size =  task_p_p + 1;
	MPI_Request req[2];
	MPI_Status st[2];
	for (int i = 0; i < p_count; ++i)
	{		
		block_i = (m_rank + i) % p_count;
		multiMatVecAdd(l_matrix , l_vec_repr , l_vec_dest , l_matrix_r , recvcounts[block_i]  , l_matrix_c , displs[block_i]);
		if(p_count > 1){		
			MPI_Isend(l_vec_repr , exc_size ,MPI_DOUBLE, (m_rank-1+p_count)%p_count ,12345,MPI_COMM_WORLD,&req[0]);
			MPI_Irecv(l_vec_recv_repr , exc_size , MPI_DOUBLE, (m_rank + 1)%p_count ,12345,MPI_COMM_WORLD,&req[1]);
			MPI_Waitall(2,req,st);
			double* keep = l_vec_repr;
			l_vec_repr = l_vec_recv_repr;
			l_vec_recv_repr = keep;
		}
	}
	if(p_count > 1)
		copyVector(l_vec_ , l_vec_repr , l_matrix_r);
}





void absolute_value(double* l_vec_ , double* dest ,  int size){
	double sum = 0;
	for (int i = 0; i < size; ++i){
		sum+=(l_vec_[i] * l_vec_[i]);
	}
	MPI_Reduce(&sum , dest , 1 , MPI_DOUBLE , MPI_SUM, 0, MPI_COMM_WORLD);
	*dest = sqrt(*dest);
}



void fillVectorM3(double* v , int size){
	for (int i = 0; i < size; ++i)
	{
		int r = rand();
		if(r % 2 == 0)
			v[i] = r % 100 - 50;
		else v[i] = 0;
	}
}

void fillVectorM3Grid(double* grid , double* b , int size){
	for (int i = 0; i < size; ++i)
	{
		grid[i] = b[i];
	}
}

void fillMatrixM3(double* m , int r , int c , int l_s_index , int N_x){

	for(int i = 0; i < r; ++i){
		for(int j = 0; j < c; ++j){
			if(j == l_s_index + i){
				m[i * c + j] = -4;
			}else if(j + 1 == l_s_index + i){
				if(i && i % N_x == 0){
					m[i * c + j] = 0;
				}else{
					m[i * c + j] = 1;
				}
			}else if(j - 1 == l_s_index + i){
				if(j && j % N_x == 0){
					m[i * c + j] = 0;
				}else{
					m[i * c + j] = 1;
				}
			}else if(j == l_s_index + i + N_x || j == l_s_index + i - N_x){
				m[i * c + j] = 1;
			}else{
				m[i * c + j] = 0;
			}	
		}
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


	int task_per_process = MATRIX_SIZE / p_count;
	int overflow_tasks = MATRIX_SIZE % p_count;
	int l_s_index = 0 , l_e_index = 0;
	l_s_index = 0;
	l_e_index = 0;


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
	int exc_size = task_per_process + 1;
	double* l_matrix = (double*)malloc(l_matrix_c * l_matrix_r * sizeof(double));
	double* l_vec_x = (double*)malloc(sizeof(double) * exc_size);
	double* l_vec_y = (double*)malloc(sizeof(double) * exc_size);
	double* l_vec_b = (double*)malloc(sizeof(double) * exc_size);
	double* l_vec_dest = (double*)malloc(sizeof(double) * exc_size);
	double* l_vec_recv = (double*)malloc(sizeof(double) * exc_size);
	double* answer = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double VEC_B_ABSOLUTE = 1;
	/*fillVector(l_vec_x , exc_size , 0);
	fillMatrixM1(l_matrix , l_matrix_r , l_matrix_c , l_s_index);
	/////////M2
	fillVectorM2(answer , MATRIX_SIZE);
	multiMatVecTest(l_matrix , answer  , l_vec_b , l_matrix_r , l_matrix_c);
	absolute_value(l_vec_b , &VEC_B_ABSOLUTE , l_matrix_r);*/

	//TempTest
	fillVectorM3(l_vec_b , l_matrix_r);
	fillMatrixM3(l_matrix , l_matrix_r , l_matrix_c , l_s_index , N_x);
	fillVectorM3Grid(l_vec_x , l_vec_b , l_matrix_r);

	/*MPI_Barrier(MPI_COMM_WORLD);
	if(1){
		showVector(l_vec_b , l_matrix_r);
		//showVector(l_vec_x , l_matrix_r);
	}
	return 0;
	MPI_Barrier(MPI_COMM_WORLD);
	showMatrix(l_matrix , l_matrix_r , l_matrix_c);
	return 0;*/

	double start_t = MPI_Wtime() , end_t = 0;
	double t_factor = 0;
	int end = 0;
	int count = 0;
	while(1){
		++count;
		fillVector(l_vec_y , exc_size , 0);
		multi(l_matrix , l_vec_x, l_vec_y , l_vec_recv ,recvcounts , displs , l_matrix_r\
			 , l_matrix_c ,  m_rank , p_count , task_per_process);

		subtractVec(l_vec_y , l_vec_b , l_vec_y , l_matrix_r);

		double a = 0;
		absolute_value(l_vec_y , &a , l_matrix_r);


		if(0 == m_rank){
			if(a / VEC_B_ABSOLUTE < 0.00001){
 				end = 1;
				printf("Acc: %lf\n" , a / VEC_B_ABSOLUTE);
			} 	
			//printf("Acc: %lf\n" , a / VEC_B_ABSOLUTE);
		}

		MPI_Bcast(&end , 1 , MPI_INT, 0 , MPI_COMM_WORLD);

		if(end){
			//MPI_Barrier(MPI_COMM_WORLD);
			//showVector(l_vec_x , l_matrix_r);
			end_t = MPI_Wtime();
			if(0 == m_rank){
				printf("Time taken: %f Iterations: %d \n" , end_t - start_t , count);
			}

			{
				free(l_matrix); free(l_vec_x);
				free(l_vec_y); free(l_vec_dest);
				free(l_vec_recv);free(l_vec_b);
				free(answer);
			}
			MPI_Finalize();
			return 0;
		}

		fillVector(l_vec_dest , exc_size , 0);
		multi(l_matrix , l_vec_y, l_vec_dest , l_vec_recv ,recvcounts , displs\
			 , l_matrix_r , l_matrix_c ,  m_rank , p_count , task_per_process);

		double l_d1 = 0 , l_d2 = 0 , d1 = 0 , d2 = 0;
		dotProduct(l_vec_y , l_vec_dest , l_matrix_r , &l_d1 , m_rank);
		dotProduct(l_vec_dest , l_vec_dest , l_matrix_r , &l_d2 , m_rank);
		MPI_Reduce(&l_d1 , &d1 , 1 , MPI_DOUBLE , MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&l_d2 , &d2 , 1 , MPI_DOUBLE , MPI_SUM, 0, MPI_COMM_WORLD);

		if(0 == m_rank){
			t_factor = d1 / d2;
		}
		MPI_Bcast(&t_factor , 1 , MPI_DOUBLE, 0 , MPI_COMM_WORLD);
		multiOnScalar(l_vec_y , t_factor , l_matrix_r);
		subtractVec(l_vec_x , l_vec_y , l_vec_x , l_matrix_r);
	}

	return 0;

}
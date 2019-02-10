#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>

//передавать -lm

#define MATRIX_SIZE 100

//имметричные матрицы
//метод минимальных невязок


/*y(n) = A * x(n) - b;


t(n) = ( y(n) dotp (A * y(n)) ) / ( (A * y(n)) dotp ( (A * y(n)))

x(n + 1) = x(n) - t(n) * y(n)

|A * x(n)| / |b| < E = 10^(-5)



1) dotp
2) | |
3) A * x

t , x , y - отправляю всем                */



/*0----| |  |
1	 | |= |
2	 | |  |
3----| |  |
4*/


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

//неверно
/*void fillMatrix(double* m , int r , int c){
	for(int i = 0; i < r; ++i){
		for(int j = 0; j < c; ++j){
			m[i * r + j] = (j == i ? 2.0 : 1.0);
		}
	}
}*/

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
	if(0 == a2){
		 printf("|b| = 0?\n");
		 exit(0);
	}

	return a1 / a2;
}



int main(int argc , char** argv){

	MPI_Init(&argc  , &argv);

	int m_rank = 0 , p_count = 0;

	MPI_Comm_size(MPI_COMM_WORLD , &p_count);
	MPI_Comm_rank(MPI_COMM_WORLD , &m_rank);


	int task_per_process = MATRIX_SIZE / p_count;
	int overflow_tasks = MATRIX_SIZE % p_count;
	//start index ,  end index
	//идексы включительно
	int l_s_index = 0 , l_e_index = 0;
	l_s_index = m_rank * task_per_process;
	l_e_index = l_s_index + task_per_process - 1;


	int* recvcounts = (int*)malloc(sizeof(int) * p_count);
	int* displs = (int*)malloc(sizeof(int) * p_count);
	int i_sum = 0;
	for (int i = 0; i < p_count; ++i)
	{
		double s = i * task_per_process;
		double e = s + task_per_process - 1;
		if(0 != overflow_tasks){
			if(i < overflow_tasks){
				s+=i;
			    e+=(i + 1);
			}else{
				s+=overflow_tasks;
				e+=overflow_tasks;
			}
		}
		displs[i] = i_sum;
		recvcounts[i] = e - s + 1;
		i_sum+=recvcounts[i];
	}


	//можно поместить в цикл выше
	if(0 != overflow_tasks){
		if(m_rank < overflow_tasks){
			l_s_index+=m_rank;
		    l_e_index+=(m_rank + 1);
		}else{
			l_s_index+=overflow_tasks;
			l_e_index+=overflow_tasks;
		}
	}

	int l_matrix_r = (l_e_index - l_s_index + 1);
	int l_matrix_c = MATRIX_SIZE;
	double* l_matrix = (double*)malloc(l_matrix_c * l_matrix_r * sizeof(double));
	double* l_vec_keepAx = (double*)malloc(sizeof(double) * l_matrix_r);


	fillMatrix(l_matrix , l_matrix_r , l_matrix_c);


	double* vec_x = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_b = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_Gathered_y = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_Gathered = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double  t_factor = 0;
	fillVector(vec_x , MATRIX_SIZE , 0);
	fillVector(vec_b , MATRIX_SIZE , MATRIX_SIZE + 1);

	//освободить память


	printf("%d -  ls: %d , le: %d\n", m_rank , l_s_index , l_e_index);
	printf("r: %d , c: %d --\n", l_matrix_r , l_matrix_c );

	//MPI_GatherV

	while(1){

		multiMatVec(l_matrix ,  vec_x , l_vec_keepAx , l_matrix_r , l_matrix_c);

		MPI_Allgatherv(l_vec_keepAx , l_matrix_r , MPI_DOUBLE,vec_Gathered ,
		recvcounts , displs ,  MPI_DOUBLE, MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD);

		//printf("End\n");

		//calculate y(n)
		subtractVec(vec_Gathered , vec_b , vec_Gathered_y , MATRIX_SIZE);

		if(0 == m_rank){
			double check = checkAccuracy(vec_Gathered_y , vec_b , MATRIX_SIZE);
			if(check < 0.00001){
				showVector(vec_x , MATRIX_SIZE);
				printf("-%lf\n", check);	
				return 0;
			}else{
				//printf("----------%lf\n", check );
			}
		}

		//calculate A * y(n)
		multiMatVec(l_matrix , vec_Gathered_y , l_vec_keepAx , l_matrix_r , l_matrix_c);

		MPI_Allgatherv(l_vec_keepAx , l_matrix_r , MPI_DOUBLE,vec_Gathered ,
		recvcounts , displs ,  MPI_DOUBLE, MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD);

		double dotpr = 0;
		dotProduct(vec_Gathered_y , vec_Gathered , MATRIX_SIZE , &dotpr);
		t_factor = dotpr;
		dotProduct(vec_Gathered , vec_Gathered , MATRIX_SIZE , &dotpr);
		t_factor/=dotpr;

		multiOnScalar(vec_Gathered_y , t_factor , MATRIX_SIZE);

		subtractVec(vec_x , vec_Gathered_y , vec_x , MATRIX_SIZE);
		//printf("Hello\n");
	}





	MPI_Finalize();
	return 0;

}
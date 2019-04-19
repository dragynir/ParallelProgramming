#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>

#define MATRIX_SIZE 10


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

void addVec(double* a , double* b , double* dest , int size){
	for(int i = 0; i < size; ++i){
		dest[i] = (a[i] + b[i]);
	}
}



void multiOnScalar(double* v , double scalar , int size){
	for (int i = 0; i < size; ++i)
	{
		v[i]*=scalar;
	}
}

void fillMatrix(double* m , int r , int c , int l_s_index){

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


void absoluteValue(double *a , int size ,  double* dest){
	double d = 0;
	dotProduct(a , a , size , &d , 1);
	*dest = sqrt(d);
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

void copyVector(double* a , double* b , int size){
	for (int i = 0; i < size; ++i)
	{
		a[i] = b[i];
	}
}

void fillVectorM2(double* v , int size){
	for (int i = 0; i < size; ++i)
	{
		v[i] = 0.7;
		//v[i] = sin((3.14 * 2 * i) / size);
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

int main(int argc , char** argv){

	MPI_Init(&argc  , &argv);

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
	fillVector(l_vec_x , exc_size , 0);
	fillMatrix(l_matrix , l_matrix_r , l_matrix_c , l_s_index);
	double* dotpResiver = (double*)malloc(sizeof(double) * p_count * 2);
	double* dotpSender = (double*)malloc(sizeof(double) * 2);
	double* vec_b = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_Gathered = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double VEC_B_ABSOLYTE = 1;

	//M1
	//fillVector(vec_b , MATRIX_SIZE , MATRIX_SIZE + 1);
	//fillVector(l_vec_b ,  exc_size , MATRIX_SIZE + 1);
	//absoluteValue(vec_b , MATRIX_SIZE , &VEC_B_ABSOLYTE);

	/////////M2
	fillMatrix(l_matrix , l_matrix_r , l_matrix_c , l_s_index);
	fillVectorM2(vec_Gathered , MATRIX_SIZE);
	multiMatVecTest(l_matrix ,  vec_Gathered , l_vec_recv , l_matrix_r , l_matrix_c);
	MPI_Allgatherv(l_vec_recv , l_matrix_r , MPI_DOUBLE,vec_b ,
	recvcounts , displs ,  MPI_DOUBLE, MPI_COMM_WORLD);
	absoluteValue(vec_b , MATRIX_SIZE , &VEC_B_ABSOLYTE);
	/////////M2

/*y(n) = A * x(n) - b;

t(n) = ( y(n) dotp (A * y(n)) ) / ( (A * y(n)) dotp ( (A * y(n)))

x(n + 1) = x(n) - t(n) * y(n)

|A * x(n) - b| / |b| < E = 10^(-5)*/


	/*printf("%d -  ls: %d , le: %d\n", m_rank , l_s_index , l_e_index);
	printf("r: %d , c: %d --\n", l_matrix_r , l_matrix_c );
	if(0 == m_rank){
		for (int i = 0; i < p_count; ++i)
		{
			printf("displs: %d recvcounts: %d\n", displs[i] , recvcounts[i]);
		}
	}*/

	MPI_Request req[2];
	MPI_Status st[2];
	double t_factor = 0;
	int count = 0;
	while(1){
		//Request?	

		fillVector(l_vec_dest , exc_size , 0);
		double* l_vec_repr = l_vec_x;
		double* l_vec_recv_repr = l_vec_recv;
		int block_i = 0;
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
			copyVector(l_vec_x , l_vec_repr , l_matrix_r);


		subToTest(l_vec_dest , vec_b , l_vec_y, l_matrix_r , l_s_index);
		//subtractVec(l_vec_dest , l_vec_b , l_vec_y , l_matrix_r);
	
		MPI_Allgatherv(l_vec_y , l_matrix_r , MPI_DOUBLE,vec_Gathered , recvcounts , displs, MPI_DOUBLE ,MPI_COMM_WORLD);

		double a = 0;
		absoluteValue(vec_Gathered , MATRIX_SIZE , &a);

		if(a / VEC_B_ABSOLYTE < 0.0001){
			MPI_Gatherv(l_vec_x , l_matrix_r , MPI_DOUBLE,vec_Gathered , recvcounts , displs, MPI_DOUBLE , 0 , MPI_COMM_WORLD);
			if(0 == m_rank){
				showVector(vec_Gathered , MATRIX_SIZE);
				printf("Acc: %lf , %d\n" , a / VEC_B_ABSOLYTE , count);
			}
			return 0;
		}

		fillVector(l_vec_dest , exc_size , 0);
		l_vec_repr = l_vec_y;
		l_vec_recv_repr = l_vec_recv;
		for (int i = 0; i < p_count; ++i)
		{
			block_i = (m_rank + i) % p_count;
			multiMatVecAdd(l_matrix , l_vec_repr , l_vec_dest , l_matrix_r , recvcounts[block_i] , l_matrix_c , displs[block_i]);
			if(p_count > 1){	
				MPI_Isend(l_vec_repr , exc_size ,MPI_DOUBLE, (m_rank-1+p_count)%p_count ,12340,MPI_COMM_WORLD,&req[0]);
				MPI_Irecv(l_vec_recv_repr , exc_size , MPI_DOUBLE, (m_rank + 1)%p_count  ,12340,MPI_COMM_WORLD,&req[1]);
				MPI_Waitall(2,req,st);
				double* keep = l_vec_repr;
				l_vec_repr = l_vec_recv_repr;
				l_vec_recv_repr = keep;
			}
		}
		if(p_count > 1)
			copyVector(l_vec_y , l_vec_repr , l_matrix_r);
		//l_vec_y = l_vec_repr;
		

		dotProduct(l_vec_y , l_vec_dest , l_matrix_r , dotpSender , m_rank);
		dotProduct(l_vec_dest , l_vec_dest , l_matrix_r ,dotpSender + 1 , m_rank);
		MPI_Gather(dotpSender, 2 , MPI_DOUBLE , dotpResiver, 2, MPI_DOUBLE, 0 , MPI_COMM_WORLD);

		if(0 == m_rank){
			double d1 = 0 , d2 = 0;
			for (int i = 0; i < p_count * 2; ++i)
			{
				if(i % 2 == 0) d1+=dotpResiver[i];
				else d2+=dotpResiver[i];
			}
			t_factor = d1 / d2;
		}

		MPI_Bcast(&t_factor , 1 , MPI_DOUBLE, 0 , MPI_COMM_WORLD);
		multiOnScalar(l_vec_y , t_factor , l_matrix_r);
		subtractVec(l_vec_x , l_vec_y , l_vec_x , l_matrix_r);
	}

	MPI_Finalize();
	return 0;

}
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<sys/time.h>
#include<time.h>
#include<omp.h>



void dotProduct(double* a , double* b , int size , double* dest){
	#pragma omp for
	for(int i = 0; i < size; ++i){
		(*dest)+=(a[i] * b[i]);
	}
}



void absoluteValue(double *a , int size ,  double* dest){
	double d = 0;
	#pragma omp for
	for(int i = 0; i < size; ++i){
		(*dest)+=(a[i] * a[i]);
	}
	*dest = sqrt(d);
}


void multiMatVec(double* A ,  double* x , double* dest , int A_rows , int A_col){
	#pragma omp for 
	for(int r = 0; r < A_rows; ++r){
		double s = 0;
		for(int c = 0; c < A_col; ++c){
			s+=(A[r * A_col + c] * x[c]); 
		}
		dest[r] = s;
	}
}

void subtractVec(double* a , double* b , double* dest , int size){
	#pragma omp for
	for(int i = 0; i < size; ++i){
		dest[i] = (a[i] - b[i]);
	}
}



void fillVector(double* v , int size , double value){
	#pragma omp for
	for (int i = 0; i < size; ++i)
	{
		v[i] = value;
	}
}


void showVector(double* v , int size){
	#pragma omp for
	for (int i = 0; i < size; ++i)
	{
		printf("%lf\n", v[i]);
	}
}


double checkAccuracy(double* Ax , double abs_b , int size){
	double a1 = 0 , a2 = 0;
	absoluteValue(Ax , size , &a1);
	if(0 == abs_b){
		 printf("|b| = 0?\n");
		 exit(0);
	}

	return a1 / abs_b;
}

void fillVectorM3(double* v , int size){
	#pragma omp  for
	for (int i = 0; i < size; ++i)
	{
		v[i] = 1;
		if(!(i%3)) v[i] = (1 - 2*(i%2))*50; 
		else v[i] = 0;
	}
}


void fillMatrixM3(double* m , int r , int c , int N_x){
	#pragma omp for
	for(int i = 0; i < r; ++i){
		for(int j = 0; j < c; ++j){
			if(j == i){
				m[i * c + j] = -4;
			}else if(j == i + N_x || j == i - N_x){
				m[i * c + j] = 1;
			}else if(j == i + 1){
				if(j && (j % N_x == 0)){
					m[i * c + j] = 0;
				}else{
					m[i * c + j] = 1;
				}
			}else if(j == i - 1){
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
	#pragma omp for
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			printf("It: %d %d , %d", i , j , (int)m[i * c + j]);
		}
		printf("\n");
	}
}



/*void calculateTFactor(double* a1 , double* b1 ,  double* a2 , double* b2 ,  int size , double* s1 , double* s2){
	#pragma omp for reduction(+:*s1,*s2)
	for(int i = 0; i < size; ++i){
		s1+=(a1[i] * b1[i]);
		s2+=(a2[i] * b2[i]);
	}
}*/

void calculateX(double* vec_x , double* vec_y ,  int size , double t_factor){
	#pragma omp for
	for (int i = 0; i < size; ++i)
	{
		vec_x[i] = vec_x[i] - vec_y[i] * t_factor;
	}
}



void fillMatrixM1(double* m , int r , int c){
	#pragma omp for
	for(int i = 0; i < r; ++i){
		for(int j = 0; j < c; ++j){
			if((j == i)){
				m[i * c + j] = 2;
			}else{
				m[i * c + j] = 1;
			}
		}
	}
}


/*gcc –fopenmp –o hello hello.c
> OMP_NUM_THREADS=4 ./hello*/

int foo(double* v , int size){
	#pragma omp for
	for(int i = 0; i < size; ++i)
	{
		printf("%lf\n", 1.0);
	}
	return 2;
}

int get(){
	return foo(0 , 1);
}

int main(int argc , char** argv){
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
	double* matrix_a = (double*)malloc(sizeof(double) * MATRIX_SIZE * MATRIX_SIZE);
	double* vec_x = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_b = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_keep = (double*)malloc(sizeof(double) * MATRIX_SIZE);
	double* vec_y = (double*)malloc(sizeof(double) * MATRIX_SIZE);	

	#pragma omp parallel
	{
		/////////TempTest
		fillVectorM3(vec_b , MATRIX_SIZE);
		fillMatrixM3(matrix_a , MATRIX_SIZE , MATRIX_SIZE , N_x);
		fillVector(vec_x , MATRIX_SIZE , 0);

		/////M1
		/*fillMatrixM1(matrix_a , MATRIX_SIZE , MATRIX_SIZE);
		fillVector(vec_x , MATRIX_SIZE, 0.0002);
		multiMatVec(matrix_a ,  vec_x , vec_b , MATRIX_SIZE , MATRIX_SIZE);
		fillVector(vec_x , MATRIX_SIZE , 0);*/
	}



	struct timeval start, end;
	gettimeofday(&start, NULL);

	
	double absolute_b = 0 , absolute_vec = 0;
	double s1 = 0 , s2 = 0;
	int count = 0;
	#pragma omp parallel
	{

		#pragma omp for reduction(+:absolute_b)
		for(int i = 0; i < MATRIX_SIZE; ++i){
			absolute_b+=(vec_b[i] * vec_b[i]);
		}
		#pragma omp single
		absolute_b = sqrt(absolute_b);



		//showVector(vec_b , MATRIX_SIZE);
		//printf("%lf\n", absolute_b);

		//showMatrix(matrix_a , MATRIX_SIZE , MATRIX_SIZE);
		//showVector(vec_x , MATRIX_SIZE);
		//showVector(vec_b , MATRIX_SIZE);

		while(1){
			//break;
			#pragma omp single //nowait
			++count;
			
			multiMatVec(matrix_a ,  vec_x , vec_keep , MATRIX_SIZE , MATRIX_SIZE);
			subtractVec(vec_keep , vec_b , vec_y , MATRIX_SIZE);
			#pragma omp single
			absolute_vec = 0;

			#pragma omp for reduction(+:absolute_vec)
			for(int i = 0; i < MATRIX_SIZE; ++i){
				absolute_vec+=(vec_y[i] * vec_y[i]);
			}
			double check = sqrt(absolute_vec) / absolute_b;
			if(check < 0.00001){
				printf("res: %lf\n", check);
				break;
			}else{
				//printf("%lf\n", check);
			}

			multiMatVec(matrix_a , vec_y , vec_keep	 , MATRIX_SIZE , MATRIX_SIZE);
			#pragma omp single
			{
				s1 = 0;
				s2 = 0;
			}

			/*#pragma omp for reduction(+:s1,s2)
			for(int i = 0; i < MATRIX_SIZE; ++i){
				s1+=(vec_y[i] * vec_keep[i]);
				s2+=(vec_keep[i] * vec_keep[i]);
			}*/
			#pragma omp for reduction(+:s1)
			for(int i = 0; i < MATRIX_SIZE; ++i){
				s1+=(vec_y[i] * vec_keep[i]);
			}

			#pragma omp for reduction(+:s2)
			for(int i = 0; i < MATRIX_SIZE; ++i){
				s2+=(vec_keep[i] * vec_keep[i]);
			}
			double t_factor = s1 / s2;
			calculateX(vec_x , vec_y , MATRIX_SIZE , t_factor);
		}

	}




	
	//clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	gettimeofday(&end, NULL);

	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
 		end.tv_usec - start.tv_usec) / 1.e6;
		//printf("Time taken: %lf sec. ",
				//end.tv_sec-start.tv_sec + 0.000000001*(end.tv_nsec-start.tv_nsec));
	printf("Time taken: %lf sec. Iterations: %d \n" ,delta ,  count);

    //showVector(vec_x , MATRIX_SIZE);
	{
	free(matrix_a);free(vec_keep);
	free(vec_x);free(vec_b);
	free(vec_y);
	}

	

}

/*/*y(n) = A * x(n) - b;


t(n) = ( y(n) dotp (A * y(n)) ) / ( (A * y(n)) dotp ( (A * y(n)))

x(n + 1) = x(n) - t(n) * y(n)

|A * y(n)| / |b| < E = 10^(-5)*/


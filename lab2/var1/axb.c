#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<sys/time.h>
#include<time.h>
#include<omp.h>



void dotProduct(double* 
	a , double* b , int size , double* dest){
	double s = 0;
	#pragma omp parallel for reduction(+:s)
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
	#pragma omp parallel for 
	for(int r = 0; r < A_rows; ++r){
		double s = 0;
		for(int c = 0; c < A_col; ++c){
			s+=(A[r * A_col + c] * x[c]); 
		}
		dest[r] = s;
	}
}


void subtractVec(double* a , double* b , double* dest , int size){
	#pragma omp parallel for
	for(int i = 0; i < size; ++i){
		dest[i] = (a[i] - b[i]);
	}
}



void fillVector(double* v , int size , double value){
	#pragma omp parallel for
	for (int i = 0; i < size; ++i)
	{
		v[i] = value;
	}
}


void showVector(double* v , int size){
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
	#pragma omp parallel for
	for (int i = 0; i < size; ++i)
	{
		v[i] = 1;
		if(!(i%3)) v[i] = (1 - 2*(i%2))*50; 
		else v[i] = 0;
	}
}


void fillMatrixM3(double* m , int r , int c , int N_x){
	#pragma omp parallel for
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
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			printf("%d ", (int)m[i * c + j]);
		}
		printf("\n");
	}
}



void calculateTFactor(double* a1 , double* b1 ,  double* a2 , double* b2 ,  int size , double* dest){
	double s1 = 0 , s2 = 0;
	#pragma omp parallel for reduction(+:s1,s2)
	for(int i = 0; i < size; ++i){
		s1+=(a1[i] * b1[i]);
		s2+=(a2[i] * b2[i]);
	}
	*dest = s1 / s2;
}

void calculateX(double* vec_x , double* vec_y ,  int size , double t_factor){
	#pragma omp parallel for
	for (int i = 0; i < size; ++i)
	{
		vec_x[i] = vec_x[i] - vec_y[i] * t_factor;
	}
}



void fillMatrixM1(double* m , int r , int c){

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


	/////////TempTest
	fillVectorM3(vec_b , MATRIX_SIZE);
	fillMatrixM3(matrix_a , MATRIX_SIZE , MATRIX_SIZE , N_x);
	fillVector(vec_x , MATRIX_SIZE , 0);

	/////M1
	/*fillMatrixM1(matrix_a , MATRIX_SIZE , MATRIX_SIZE);
	fillVector(vec_x , MATRIX_SIZE, 0.0001);
	multiMatVec(matrix_a ,  vec_x , vec_b , MATRIX_SIZE , MATRIX_SIZE);
	fillVector(vec_x , MATRIX_SIZE , 0);*/


	/*struct timespec start, end;
 	clock_gettime(CLOCK_MONOTONIC_RAW, &start);*/



 	struct timeval start, end;
	gettimeofday(&start, NULL);



	double absolute_b = 0;
	absoluteValue(vec_b	 , MATRIX_SIZE , &absolute_b);	
	int count = 0;
	while(1){
		++count;
		multiMatVec(matrix_a ,  vec_x , vec_keep , MATRIX_SIZE , MATRIX_SIZE);
		//calculate y(n)
		subtractVec(vec_keep , vec_b , vec_y , MATRIX_SIZE);

		double check = 1;
		check = checkAccuracy(vec_y	 , absolute_b , MATRIX_SIZE);
		if(check < 0.00001){
			printf("res: %lf\n", check);

		
			gettimeofday(&end, NULL);

			double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + 
         		end.tv_usec - start.tv_usec) / 1.e6;
 			
			printf("Time taken: %lf sec. Iterations: %d \n" ,delta ,  count);

		    
			{
			free(matrix_a);free(vec_keep);
			free(vec_x);free(vec_b);
			free(vec_y);
			}
			return 0;
		}
		
		//calculate A * y(n)
		multiMatVec(matrix_a , vec_y , vec_keep	 , MATRIX_SIZE , MATRIX_SIZE);
		double  t_factor = 0;
		calculateTFactor(vec_y	 , vec_keep	 ,  vec_keep , vec_keep	,  MATRIX_SIZE	, &t_factor);
		calculateX(vec_x , vec_y , MATRIX_SIZE , t_factor);
	}

	return 0;

}


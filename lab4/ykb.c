#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#define E 10e-8
#define A 10e+5

//512 384 256 : 36.6 sec
//Косяк при 3 потока 4 4 4 и их четных братьев

//Сигфолт при MPI_Finalize()

double fi(double x, double y, double z){
	return x*x + y*y + z*z;
}

double ro(double x, double y, double z){
	return	 6 - A * fi(x,y,z);
}



void defineSlices(int* displs, int* counts, int per_proc,
	 int over, int p_count,  int m_rank, int* start, int* end)
{
		int i_sum = 0;
		for (int i = 0; i < p_count; ++i)
		{
			double s = i * per_proc + (i < over ? i : over);
			double e = s + per_proc - 1 + (i < over ? 1 : 0);

			displs[i] = i_sum;
			counts[i] = e - s + 1;
			i_sum+=counts[i];

			if(i == m_rank){
				*start = s;
				*end = e;
			}
		}
}


//Функция заполняет слой
// Каждому узлу присваивает конкретную точку
//Иницилизирует начальные значения
void fillSlice(double* l_slice, double* model_space, double* space_displ, 
 int* node_counts, double slice_start_z, double slice_end_z, int p_count)

{
	//заполнение идет с верней грани слоя
	double curr_z_start = model_space[5] -  slice_start_z * space_displ[2];
	/*printf("%lf\n", curr_z_start );
	return;*/

	int slice_z = slice_end_z - slice_start_z + 1;
	int slice_x = node_counts[0];
	int slice_y = node_counts[1];
	//printf("%d %d %d\n", slice_x , slice_y , slice_z);
	double value = 0;
	for (int i = 0; i < slice_z; ++i)
	{
		double curr_z =curr_z_start - space_displ[2] * i;
		/*заполнение идет с нижнего левого угла
				^y
				|
		________|________
		|		|		|
		|		|		|
		|_______|_______|____> x
		|		|		|
		|		|		|
	   *|_______|_______|                 */


		double curr_x = 0;
		double curr_y = 0;
		//цикл по x
		for (int j = 0; j < slice_x; ++j)
		{
			curr_x = space_displ[0] * j + model_space[0];
			//цикл по y
			for (int k = 0; k < slice_y; ++k)
			{
				curr_y = space_displ[1] * k + model_space[2];
	

				//printf("(%lf , %lf , %lf) ", curr_x, curr_y, curr_z);

				if(slice_start_z == 0 && i == 0){
					value = fi(curr_x , curr_y , curr_z);
					//printf("+\n");
				}else if(slice_end_z == node_counts[2] - 1 && i == slice_z - 1){
					value = fi(curr_x , curr_y , curr_z);
					//printf("+\n");
				}
				else if(k == 0 || k == slice_y - 1 || j == 0 || j == slice_x - 1){
					value = fi(curr_x , curr_y , curr_z);
					//printf("+\n");
				}else{
					value = 0;//значение по умолчанию
				}
				l_slice[i * slice_x * slice_y + j * slice_y + k] = value;
			}
			//printf("\n");
		}
		//printf("\n\n");
	}
}



double getMaxDiff(double* l_slice, double* keep_l_slice, int slice_x, int slice_y, int slice_z){
	double max = 0;
	for (int i = 0; i < slice_z; ++i)
	{
		for (int j = 0; j < slice_x; ++j)
		{
			for (int k = 0; k < slice_y; ++k)
			{
				int index = i * slice_x * slice_y + j * slice_y + k;
				double n = l_slice[index]; 
				double n_1 = keep_l_slice[index];
				double diff =  fabs(n - n_1);
				if(diff > max){
					max = diff;
				}
			}	
		}
	}
	return max;
}


void showSlice(double* l_slice, int slice_x, int slice_y , int slice_z){
	for (int i = 0; i < slice_z; ++i)
	{
		for (int j = 0; j < slice_x; ++j)
		{
			for (int k = 0; k < slice_y; ++k)
			{
				double* n = l_slice + i * slice_x * slice_y + j * slice_y + k; 
				printf("%lf ", *n);
			}
			printf("\n");	
		}
		printf("\n\n");	
	}
}





void doIteration(double* l_slice, double* keep_l_slice, int z_center, int slice_x, int slice_y,
	double hx, double hy, double hz, double MULTI_KF, double* model_space){


	/*Тут нет проверки на выход за границы зоны->
	сигфолты велком.Непонятно как использовать формулу на границе*/
	//printf("%lf\n", MULTI_KF);
	//printf("%lf , %lf , %lf\n", hx , hy , hz );
	//int left_edge = 0 , right_edge = 

	double z = model_space[5] -  z_center * hz;

	for (int i = 0; i < slice_x; ++i)
		{
			for (int j = 0; j < slice_y; ++j)
			{
				int dest = z_center *  slice_x * slice_y + i * slice_y + j;
				if(i == 1 || j == 1 || j == slice_y - 1 || i == slice_x){
					keep_l_slice[dest] = l_slice[dest];
					continue;
				}
				double fi_z = (l_slice[dest + slice_x * slice_y] +
				 		l_slice[dest - slice_x * slice_y]) / (hz * hz);
				double fi_x = (l_slice[dest + slice_y] + l_slice[dest - slice_y]) /  (hx * hx);

				double fi_y = (l_slice[dest + 1] + l_slice[dest - 1]) / (hy * hy);
				double x = hx * i + model_space[0];
				double y = hy * j + model_space[2];
				keep_l_slice[dest] = MULTI_KF * ((fi_x + fi_y + fi_z) - ro(x,y,z));
			}
		}
}




int main(int argc , char** argv){
	MPI_Init(&argc  , &argv);
	if(argc != 4){
		printf("Wrong args\n");
		return 0;
	}
	double model_space[] = {-1 , 1 , -1 , 1 , -1 , 1};
	double space_distance = 2;

	int Nx = atoi(argv[1]);
	int Ny = atoi(argv[2]);
	int Nz = atoi(argv[3]);

	int node_counts[] = {Nx , Ny , Nz};
	double space_displ[] = {space_distance / (Nx - 1) , space_distance / (Ny - 1) , space_distance / (Nz - 1)};
	double hx = space_displ[0];
	double hy = space_displ[1];
	double hz = space_displ[2];
	double MULTI_KF = 1 / ( (2 / (hx * hx)  +  2 / (hy * hy)  +  2 / (hz * hz)) + A);
	

	int p_count , m_rank;
	MPI_Comm_size(MPI_COMM_WORLD , &p_count);
	MPI_Comm_rank(MPI_COMM_WORLD , &m_rank);

	int* counts = (int*)calloc(p_count , sizeof(int));
	int* displs = (int*)calloc(p_count , sizeof(int));
	int per_proc = Nz / p_count;
	int over = Nz % p_count;

	int slice_start_z , slice_end_z , slice_z;
	defineSlices(displs, counts, per_proc, over, p_count, m_rank, &slice_start_z, &slice_end_z);

	



	//дополнительные слои обмена
	if(m_rank != 0){
		slice_start_z-=1;
	}
	if(m_rank != p_count - 1){
		slice_end_z+=1;
	}



	slice_z = slice_end_z  - slice_start_z + 1;

	//printf("%d\n", slice_z );
	//return 0;

	double* l_slice = (double*)malloc(Nx * Ny * slice_z * sizeof(double));
	double* keep_l_slice = (double*)malloc(Nx * Ny * slice_z * sizeof(double));
	



	//задаем значение каждому узлу
	fillSlice(l_slice, model_space, space_displ, node_counts, slice_start_z, slice_end_z ,p_count);

	//сужаем область моделирования
	int slice_x = node_counts[0];
	int slice_y = node_counts[1];

	/*if(m_rank ==0)
	showSlice(l_slice, slice_x, slice_y , slice_z);
	return 0;*/

	int count = 0;

	/*printf("rank: %d , z_s: %d , z_e: %d ,  slice_z: %d,\n",m_rank
	 , slice_start_z , slice_end_z , slice_z);*/
	//return 0;
	double* fi_diff = (double*)malloc(sizeof(double) * (p_count));

	double start_t = MPI_Wtime() , end_t = 0;
	while(1){
		//printf("LOL1\n");
		++count;
		
		//TODO вычисление боковых значений
		/*if(count == 2){
			if(m_rank ==0)
				showSlice(l_slice, slice_x, slice_y , slice_z);
			break;
		}*/
		//считаем боковые значения

		//если каждый слой по 1 много пересчета добавить проверки
		doIteration(l_slice, keep_l_slice, 1, slice_x, slice_y , hx , hy , hz, MULTI_KF, model_space);
		doIteration(l_slice, keep_l_slice,  slice_z - 2, slice_x, slice_y , hx , hy , hz, MULTI_KF, model_space);
		//if(m_rank ==0)
		//	showSlice(keep_l_slice, slice_x, slice_y , slice_z);
		//return 0;



		//------------------OK


		//TODO асинхронная передача


		//Передача не MPI_DOUBLE ,  double



		if(slice_start_z == 0){
			memcpy(keep_l_slice , l_slice , slice_x * slice_y * sizeof(double));
		}else if(slice_end_z == node_counts[2] - 1){
			memcpy(keep_l_slice + (slice_z - 1) *  slice_x * slice_y
			  , l_slice + (slice_z - 1) *  slice_x * slice_y, slice_x * slice_y *sizeof(double));
		}

		//if(m_rank ==0)
		//showSlice(keep_l_slice, slice_x, slice_y , slice_z);
		//return 0;
		
		MPI_Request requests[2];
		if(p_count == 1){
			memcpy(keep_l_slice , l_slice , slice_x * slice_y * sizeof(double));
			memcpy(keep_l_slice + (slice_z - 1) *  slice_x * slice_y
			  , l_slice + (slice_z - 1) *  slice_x * slice_y, slice_x * slice_y *sizeof(double));
		}else{
			//printf("Done\n");
				///////////
			//если не последний ранг
			//Передача вниз
			if(m_rank != p_count - 1){
				//printf("Done\n");
			MPI_Isend(keep_l_slice + (slice_z - 2) *  slice_x * slice_y, slice_x * slice_y, MPI_DOUBLE,
			  m_rank + 1, 1324, MPI_COMM_WORLD, &requests[0]);
			}

			//если не первый процесс
			if(m_rank != 0){
				//printf("Done\n");
			MPI_Irecv(keep_l_slice, slice_x * slice_y, MPI_DOUBLE,
			  m_rank - 1, 1324, MPI_COMM_WORLD, &requests[0]);
			}


		
	///////////

			//если не первый процесс
			if(m_rank != 0){
			MPI_Isend(keep_l_slice + slice_x * slice_y, slice_x * slice_y, MPI_DOUBLE,
			  m_rank - 1, 1325, MPI_COMM_WORLD, &requests[1]);
			}

			//если не последний процесс
			if(m_rank != p_count - 1){
			MPI_Irecv(keep_l_slice + (slice_z - 1) *  slice_x * slice_y, slice_x * slice_y, MPI_DOUBLE,
			  m_rank + 1, 1325, MPI_COMM_WORLD, &requests[1]);
			}

	///////////
		}
	





		//showSlice(keep_l_slice, slice_x, slice_y , slice_z);
		//printf("-----------------------\n");
		//TODO выполнение локальных вычислений

		for (int i = 2 ; i < (slice_z - 1); ++i)
		{	
			doIteration(l_slice, keep_l_slice, i, slice_x, slice_y , hx , hy , hz, MULTI_KF, model_space);
		}

		


		if(p_count != 1){
			MPI_Status statuses[2];
			MPI_Waitall(2, requests, statuses);
		}
		


		//if(1 == m_rank)
		//	showSlice(keep_l_slice, slice_x, slice_y , slice_z);
		//return 0;

		//showSlice(keep_l_slice, slice_x, slice_y , slice_z);
		//return 0;



		double diff = getMaxDiff(l_slice, keep_l_slice, slice_x, slice_y, slice_z);
		int end = 0;


		MPI_Gather(&diff, 1, MPI_DOUBLE, fi_diff, p_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		if(m_rank == 0){
			for (int i = 0; i < p_count; ++i)
			{
				if(fi_diff[i] > diff){
					diff = fi_diff[i];
				}
			}
			if(diff < E){
				end = 1;
			}
		}


		//printf("LOL2\n");
		MPI_Bcast(&end , 1 , MPI_INT, 0 , MPI_COMM_WORLD);
		


		if(m_rank == 0){
			printf("%lf --- %d\n", diff , count);

		}
		if(end){
			if(m_rank == 0){
				printf("%lf\n", diff);
				/*showSlice(l_slice, slice_x, slice_y , slice_z);
				printf("-------------------\n");
				showSlice(keep_l_slice, slice_x, slice_y , slice_z);*/
			}
			end_t = MPI_Wtime();
			if(0 == m_rank){
				printf("Time taken: %f Iterations: %d \n" , end_t - start_t , count);
			}
			break;
		}



		//return 0;

		//теперь в keep_l_slice все ок

		double* ptr = keep_l_slice;
		keep_l_slice = l_slice;
		l_slice = ptr;

	}

	free(l_slice);
	free(keep_l_slice);
	MPI_Finalize();

		//TODO проверка точности между вновь вычесленными занчениями и реальной функцией


		//TODO выполнение локальных вычислений


		//TODO асинхронная передача
	return 0;
}


#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>



//Request ЧТОТОНЕТАК

#define E 10e-8
#define A 10e+5

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


void fillSlice(double* l_slice, double* model_space, double* space_displ, 
 int* node_counts, double slice_start_z, double slice_end_z, int p_count)

{
	//заполнение идет с верней грани слоя
	double curr_z_start = model_space[5] -  slice_start_z * space_displ[2];


	int slice_z = slice_end_z - slice_start_z + 1;
	int slice_x = node_counts[0];
	int slice_y = node_counts[1];
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
		for (int j = 0; j < slice_x; ++j)
		{
			curr_x = space_displ[0] * j + model_space[0];
			for (int k = 0; k < slice_y; ++k)
			{
				curr_y = space_displ[1] * k + model_space[2];
	
				if(slice_start_z == 0 && i == 0){
					value = fi(curr_x , curr_y , curr_z);
				}else if(slice_end_z == node_counts[2] - 1 && i == slice_z - 1){
					value = fi(curr_x , curr_y , curr_z);
				}
				else if(k == 0 || k == slice_y - 1 || j == 0 || j == slice_x - 1){
					value = fi(curr_x , curr_y , curr_z);
				}else{
					value = 0;
				}
				l_slice[i * slice_x * slice_y + j * slice_y + k] = value;
			}
		}
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

	

	double* l_slice = (double*)malloc(Nx * Ny * slice_z * sizeof(double));
	double* keep_l_slice = (double*)malloc(Nx * Ny * slice_z * sizeof(double));
	
	fillSlice(l_slice, model_space, space_displ, node_counts, slice_start_z, slice_end_z ,p_count);

	int slice_x = node_counts[0];
	int slice_y = node_counts[1];

	int count = 0;
	
	double start_t = MPI_Wtime() , end_t = 0;
	while(1){
		++count;
	


		doIteration(l_slice, keep_l_slice, 1, slice_x, slice_y , hx , hy , hz, MULTI_KF, model_space);
		doIteration(l_slice, keep_l_slice,  slice_z - 2, slice_x, slice_y , hx , hy , hz, MULTI_KF, model_space);
		


		if(slice_start_z == 0){
			memcpy(keep_l_slice , l_slice , slice_x * slice_y * sizeof(double));
		}else if(slice_end_z == node_counts[2] - 1){
			memcpy(keep_l_slice + (slice_z - 1) *  slice_x * slice_y
			  , l_slice + (slice_z - 1) *  slice_x * slice_y, slice_x * slice_y *sizeof(double));
		}


		MPI_Request requests[4];
		MPI_Status statuses[4];
		
		
		if(p_count == 1){
			memcpy(keep_l_slice , l_slice , slice_x * slice_y * sizeof(double));
			memcpy(keep_l_slice + (slice_z - 1) *  slice_x * slice_y
			  , l_slice + (slice_z - 1) *  slice_x * slice_y, slice_x * slice_y *sizeof(double));
		}else{

			//Передача  вверх
			if(m_rank != 0){
				MPI_Isend(keep_l_slice + slice_x * slice_y, slice_x * slice_y, MPI_DOUBLE,
			  m_rank - 1, 1325, MPI_COMM_WORLD, &requests[1]);
			}

			if(m_rank != p_count - 1){
				MPI_Irecv(keep_l_slice + (slice_z - 1) *  slice_x * slice_y, slice_x * slice_y, MPI_DOUBLE,
					m_rank + 1, 1325, MPI_COMM_WORLD, &requests[0]);
			}

			//Передача вниз
			if(m_rank != p_count - 1){
				MPI_Isend(keep_l_slice + (slice_z - 2) *  slice_x * slice_y, slice_x * slice_y, MPI_DOUBLE,
			  		m_rank + 1, 1324, MPI_COMM_WORLD, &requests[3]);
			}

			if(m_rank != 0){
				MPI_Irecv(keep_l_slice, slice_x * slice_y, MPI_DOUBLE,
			 		 m_rank - 1, 1324, MPI_COMM_WORLD, &requests[2]);
			}
		}

		
		//TODO выполнение локальных вычислений

		for (int i = 2 ; i < (slice_z - 1); ++i)
		{	
			doIteration(l_slice, keep_l_slice, i, slice_x, slice_y , hx , hy , hz, MULTI_KF, model_space);
		}

			

		if(p_count != 1){
			if(m_rank != p_count - 1){
				MPI_Wait(requests , statuses);
			}
			if(m_rank != 0){
				MPI_Wait(requests + 1 , statuses + 1);
			}
			if(m_rank != 0){
				MPI_Wait(requests + 2 , statuses + 2);
			}
			if(m_rank != p_count - 1){
				MPI_Wait(requests + 3 , statuses + 3);
			}
		}


		double min_diff , diff = getMaxDiff(l_slice, keep_l_slice, slice_x, slice_y, slice_z);
		int end = 0;

		MPI_Allreduce(&diff, &min_diff, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		if(min_diff < E){
			end = 1;
		}

		if(end){
			end_t = MPI_Wtime();
			if(0 == m_rank){
				printf("Diff: %lf\n", diff);
				printf("Proc: %d , Time taken: %f Iterations: %d \n" , p_count, end_t - start_t , count);
			}
			break;
		}

		double* ptr = keep_l_slice;
		keep_l_slice = l_slice;
		l_slice = ptr;
	}

	free(l_slice);
	free(keep_l_slice);
	free(counts);
	free(displs);


	MPI_Finalize();
	return 0;
}




#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>



void fillMatrixA(double* m , int col , int row){
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			m[i * col + j] = i * col + j;
		}
	}
}



void fillMatrixB(double* m , int col , int row){
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			m[i * col + j] = i * col + j;
		}
		
	}
}




void transpose(double* m , int col , int row){
	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			double keep = m[i * col + j];
			m[i * col + j] = m[j * col + i];
			m[j * col + i] = keep;
		}
		
	}
}



void defineVarSets(int* displs , int* sendcounts , int per_proc ,
	 int over , int rank_count ,  int coord_rank , int* start , int* end)
	{

		int i_sum = 0;
		for (int i = 0; i < rank_count; ++i)
		{
			double s = i * per_proc + (i < over ? i : over);
			double e = s + per_proc - 1 + (i < over ? 1 : 0);

			displs[i] = i_sum;
			sendcounts[i] = e - s + 1;
			i_sum+=sendcounts[i];

			if(i == coord_rank){
				*start = s;
				*end = e;
			}
		}
	}


	//прочекать
	//матрицы расширены -> умножение немного другое должно быть,
	//без учета ненужных элементов
	void multiMatrices(double* a , double* b , double* dest ,  int a_col , int a_row , int b_col){

		for (int a_r = 0; a_r < a_row; ++a_r)
		{
			for (int b_c = 0; b_c < b_col; ++b_c)
			{
				double sum = 0;
				for (int a_c = 0; a_c < a_col; ++a_c)
				{
					sum+=(a[a_r * a_col + a_c] * b[b_c * a_col + a_c]);
				}
				dest[a_r * b_col + b_c] = sum;
			}
		}
	}


int main(int argc , char** argv){
	MPI_Init(&argc  , &argv);


	if(argc != 4){
		printf("NOT e args\n");
		return 0;
	}

	const int MATRIX_A_ROW = atoi(argv[1]);
	const int MATRIX_A_COL = atoi(argv[2]);
	const int MATRIX_B_ROW = MATRIX_A_COL;
	const int MATRIX_B_COL = atoi(argv[3]);


	int dims[2]={0,0},periods[2]={0,0},coords[2],reorder=1;
	int size,rank,sizey,sizex,ranky,rankx;
	int prevy,prevx,nexty,nextx;


	MPI_Comm comm2d;
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	// определение размеров решетки: dims
	MPI_Dims_create(size,2,dims);
	sizey = dims[0]; sizex = dims[1];

	// создание коммуникатора: comm2d
	MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,reorder,&comm2d);



	// получение своих координат в двумерной решетке: coords
	MPI_Cart_get(comm2d,2,dims,periods,coords);
	ranky=coords[0]; rankx=coords[1];

	

	int *sendcounts_A = (int*)malloc(sizeof(int) * sizey);
	int *displs_A = (int*)malloc(sizeof(int) * sizey);
	int *sendcounts_B = (int*)malloc(sizeof(int) * sizex);
	int *displs_B = (int*)malloc(sizeof(int) * sizex);
	int *recvcounts	= (int*)malloc(sizeof(int) * size);
	int *displs = (int*)malloc(sizeof(int) * size);
	
	

	int splitA_per_proc = MATRIX_A_ROW / sizey; 
	int over_splitA_per_proc = MATRIX_A_ROW % sizey;
	int extended_sub_r = splitA_per_proc;
	if(over_splitA_per_proc != 0){
		extended_sub_r++;
	}
	
	int r_start_a = 0 , r_end_a = 0;

	defineVarSets(displs_A , sendcounts_A , splitA_per_proc ,
	 over_splitA_per_proc , sizey ,  ranky , &r_start_a , &r_end_a);

	int sub_a_col = MATRIX_A_COL , sub_a_row = r_end_a -  r_start_a + 1;
	double* sub_a = (double*)malloc(sizeof(double) * sub_a_col * extended_sub_r);


	int splitB_per_proc = MATRIX_B_COL / sizex; 
	int over_splitB_per_proc = MATRIX_B_COL % sizex;
	int extended_sub_c = splitB_per_proc;
	if(over_splitB_per_proc != 0){
		extended_sub_c++;
	}

	int c_start_b = 0 , c_end_b = 0;

	defineVarSets(displs_B , sendcounts_B , splitB_per_proc ,
	 over_splitB_per_proc , sizex ,  rankx , &c_start_b , &c_end_b);

	int sub_b_col = c_end_b - c_start_b + 1 , sub_b_row = MATRIX_B_ROW;
	double* sub_b = (double*)malloc(sizeof(double) * extended_sub_c * sub_b_row);




	//матрица для результата локального вычисления
	double* sub_dest = (double*)malloc(sizeof(double) * sub_a_row * sub_b_col);


	int remain_dims[2] = {0 , 1};
	MPI_Comm row_comm , col_comm;
	MPI_Cart_sub(comm2d, remain_dims, &row_comm);
	remain_dims[0] = 1;remain_dims[1] = 0;
	MPI_Cart_sub(comm2d, remain_dims, &col_comm);




	double* matrix_a = 0 , *matrix_b = 0;
	double* matrix_res = 0;

	if(rankx == 0 && ranky == 0){
		matrix_a = (double*)malloc(sizeof(double) * MATRIX_A_ROW * MATRIX_A_COL);
		matrix_b = (double*)malloc(sizeof(double) * MATRIX_B_ROW * MATRIX_B_COL);
		//matrix_res = (double*)calloc(MATRIX_A_ROW * MATRIX_B_COL , sizeof(double));
		matrix_res = (double*)malloc(MATRIX_A_ROW * MATRIX_B_COL * sizeof(double));
		fillMatrixA(matrix_a , MATRIX_A_COL ,  MATRIX_A_ROW);
		fillMatrixB(matrix_b , MATRIX_B_COL ,  MATRIX_B_ROW);
	}


	//создание типа для передачи A
	MPI_Datatype type_cont_rowA;
	MPI_Type_contiguous(MATRIX_A_COL , MPI_DOUBLE, &type_cont_rowA);
	MPI_Type_commit(&type_cont_rowA);


	//только первый столбец => rankx = 0
	if(rankx == 0){
		//Рассылка матрицы A по первому столбцу
		MPI_Scatterv(matrix_a , sendcounts_A , displs_A,
				type_cont_rowA , sub_a , sub_a_col * extended_sub_r,
				MPI_DOUBLE , 0 , col_comm);
	}



	MPI_Datatype type_vectorB;
	MPI_Aint stride = MATRIX_B_COL;
	//int stride = MATRIX_B_COL;
	MPI_Type_vector( MATRIX_B_ROW, 1 , stride, MPI_DOUBLE , &type_vectorB);

	MPI_Type_create_resized(type_vectorB,
                          0 ,
                          sizeof(double) ,
                          &type_vectorB);

	MPI_Type_commit(&type_vectorB);


	

	if(ranky == 0){
		MPI_Scatterv(matrix_b, sendcounts_B, displs,
				type_vectorB , sub_b , extended_sub_c * sub_b_row ,
	 			MPI_DOUBLE , 0 , row_comm);	
	}




	//передача матрицы sub_a
	MPI_Bcast(sub_a, sub_a_col * extended_sub_r , MPI_DOUBLE, 0 , row_comm);


	//передача матрицы sub_b
	MPI_Bcast(sub_b, extended_sub_c * sub_b_row , MPI_DOUBLE, 0, col_comm);

	//вычисление
	multiMatrices(sub_a , sub_b , sub_dest ,  sub_a_col , sub_a_row , sub_b_col);


	MPI_Datatype type_send;
	stride = MATRIX_B_COL;//!!!!!!!
	MPI_Type_vector(extended_sub_r, extended_sub_c , stride, MPI_DOUBLE , &type_send);




	MPI_Type_create_resized(type_send,
                          0 ,
                          sizeof(double) * extended_sub_c ,
                          &type_send);

	MPI_Type_commit(&type_send);

	for (int i = 0; i < size; ++i)
	{
		recvcounts[i] = 1;
	}

	for (int i = 0; i < sizey; ++i)
	{
		for (int j = 0; j < sizex; ++j)
		{
			displs[i * sizex + j] = sizex *  i * extended_sub_r + j;
		}
	} 


	MPI_Gatherv(sub_dest , sub_a_row * sub_b_col , MPI_DOUBLE	,
			matrix_res , recvcounts, displs, type_send	,
			0, comm2d);


	{
		free(sendcounts_A);
		free(displs_A);
		free(sendcounts_B);
		free(displs_B);
		free(recvcounts);
		free(displs);
	}

	if(rankx == 0 && ranky == 0){
		free(matrix_a);
		free(matrix_b);
		free(matrix_res);
	}


	
	MPI_Finalize();

	return 0;    

}

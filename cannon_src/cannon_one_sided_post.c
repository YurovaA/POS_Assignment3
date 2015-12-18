#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main (int argc, char **argv) {
	FILE *fp;
	double **A = NULL, **B = NULL, **C = NULL, *A_array = NULL, *B_array = NULL, *C_array = NULL;
	double *A_local_block = NULL, *B_local_block = NULL, *C_local_block = NULL, *A_local_block_new = NULL, *B_local_block_new = NULL;
	int A_rows, A_columns, A_local_block_rows, A_local_block_columns, A_local_block_size;
	int B_rows, B_columns, B_local_block_rows, B_local_block_columns, B_local_block_size;
	int rank, size, sqrt_size, matrices_a_b_dimensions[4];
	MPI_Comm cartesian_grid_communicator, row_communicator, column_communicator;
	MPI_Status status;
	MPI_Win win_A;
	MPI_Win win_B;

	// used to manage the cartesian grid
	int dimensions[2], periods[2], coordinates[2], remain_dims[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* For square mesh */
	sqrt_size = (int)sqrt((double) size);             
	if(sqrt_size * sqrt_size != size){
		if( rank == 0 ) perror("need to run mpiexec with a perfect square number of processes\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	// create a 2D cartesian grid 
	dimensions[0] = dimensions[1] = sqrt_size;
	periods[0] = periods[1] = 1;    
	MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &cartesian_grid_communicator);
	MPI_Cart_coords(cartesian_grid_communicator, rank, 2, coordinates);

	// create a row communicator
	remain_dims[0] = 0;            
	remain_dims[1] = 1; 
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &row_communicator);

	// create a column communicator
	remain_dims[0] = 1;
	remain_dims[1] = 0;
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &column_communicator);

	// getting matrices from files at rank 0 only
	// example: mpiexec -n 64 ./cannon matrix1 matrix2 [test]
	if (rank == 0){
		int row, column;
		// Reading A file
		if ((fp = fopen (argv[1], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[0], &matrices_a_b_dimensions[1]);
			
			// Read Matrix into A
			A = (double **) malloc (matrices_a_b_dimensions[0] * sizeof(double *));
			for (row = 0; row < matrices_a_b_dimensions[0]; row++){
				A[row] = (double *) malloc(matrices_a_b_dimensions[1] * sizeof(double));
				for (column = 0; column < matrices_a_b_dimensions[1]; column++)
					fscanf(fp, "%lf", &A[row][column]);
			}
			fclose(fp);
		} else {
			if(rank == 0) fprintf(stderr, "error opening file for matrix A (%s)\n", argv[1]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		
		// Reading B file
		// We may change rows and cols for better multiplication
		if((fp = fopen (argv[2], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[2], &matrices_a_b_dimensions[3]);
			B = (double **) malloc (matrices_a_b_dimensions[2] * sizeof(double *));
			for(row = 0; row < matrices_a_b_dimensions[2]; row++){
				B[row] = (double *) malloc(matrices_a_b_dimensions[3] * sizeof(double *));
				for(column = 0; column < matrices_a_b_dimensions[3]; column++)
					fscanf(fp, "%lf", &B[row][column]);
			}
			fclose(fp);
		} else {
			if(rank == 0) fprintf(stderr, "error opening file for matrix B (%s)\n", argv[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// need to check that the multiplication is possible given dimensions 
		// matrices_a_b_dimensions[0] = row size of A
		// matrices_a_b_dimensions[1] = column size of A
		// matrices_a_b_dimensions[2] = row size of B
		// matrices_a_b_dimensions[3] = column size of B
		if(matrices_a_b_dimensions[1] != matrices_a_b_dimensions[2]){
			if(rank == 0) fprintf(stderr, "A's column size (%d) must match B's row size (%d)\n", 
					matrices_a_b_dimensions[1], matrices_a_b_dimensions[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// this implementation is limited to cases where the matrices can be partitioned perfectly
		if( matrices_a_b_dimensions[0] % sqrt_size != 0 
				|| matrices_a_b_dimensions[1] % sqrt_size != 0 
				|| matrices_a_b_dimensions[2] % sqrt_size != 0 
				|| matrices_a_b_dimensions[3] % sqrt_size != 0 ){
			if(rank == 0) fprintf(stderr, "cannot distribute work evenly among processe\n"
					"all dimensions (A: r:%d c:%d; B: r:%d c:%d) need to be divisible by %d\n",
					matrices_a_b_dimensions[0],matrices_a_b_dimensions[1],
					matrices_a_b_dimensions[2],matrices_a_b_dimensions[3], sqrt_size );
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}

	// send dimensions to all peers
	// Cannot change Send to ISend
	if(rank == 0) {
		int i;
		for(i = 1; i < size; i++){
			MPI_Send(matrices_a_b_dimensions, 4, MPI_INT, i, 0, cartesian_grid_communicator);
		}
	} else {
		MPI_Recv(matrices_a_b_dimensions, 4, MPI_INT, 0, 0, cartesian_grid_communicator, &status);
	}

	A_rows = matrices_a_b_dimensions[0];
	A_columns = matrices_a_b_dimensions[1];
	
	B_rows = matrices_a_b_dimensions[2];
	B_columns = matrices_a_b_dimensions[3];

	// local metadata for A
	// Number of elements in each direction (local block)
	A_local_block_rows = A_rows / sqrt_size;
	A_local_block_columns = A_columns / sqrt_size;
	A_local_block_size = A_local_block_rows * A_local_block_columns;
	double* a_block = (double *) malloc (2* A_local_block_size * sizeof(double)); 
	A_local_block = a_block;
	A_local_block_new = a_block + A_local_block_size;

	// local metadata for B
	B_local_block_rows = B_rows / sqrt_size;
	B_local_block_columns = B_columns / sqrt_size;
	B_local_block_size = B_local_block_rows * B_local_block_columns;
	double* b_block = (double *) malloc (2 * B_local_block_size * sizeof(double));
	B_local_block = b_block;
	B_local_block_new = b_block + B_local_block_size;

	// We may transpose here
	
	// local metadata for C
	C_local_block = (double *) malloc (A_local_block_rows * B_local_block_columns * sizeof(double));
	// C needs to be initialized at 0 (accumulates partial dot-products)
	int i;
	for(i=0; i < A_local_block_rows * B_local_block_columns; i++){
		C_local_block[i] = 0;
	}

	// full arrays only needed at root
	if(rank == 0){
		A_array = (double *) malloc(sizeof(double) * A_rows * A_columns);
		B_array = (double *) malloc(sizeof(double) * B_rows * B_columns);
		C_array = (double *) malloc(sizeof(double) * A_rows * B_columns);
		// generate the 1D arrays of the matrices at root
		int row, column, i, j;
		// We are translating from 2d array into 1d array
		for (i = 0; i < sqrt_size; i++){
			for (j = 0; j < sqrt_size; j++){
				// Block i,j
				for (row = 0; row < A_local_block_rows; row++){
					for (column = 0; column < A_local_block_columns; column++){
						// Element inside of the block
						A_array[((i * sqrt_size + j) * A_local_block_size) + (row * A_local_block_columns) + column] 
							= A[i * A_local_block_rows + row][j * A_local_block_columns + column];
					}
				}
				for (row = 0; row < B_local_block_rows; row++){
					for (column = 0; column < B_local_block_columns; column++){
						B_array[((i * sqrt_size + j) * B_local_block_size) + (row * B_local_block_columns) + column] 
							= B[i * B_local_block_rows + row][j * B_local_block_columns + column];
					}
				}
			}
		}
		// allocate output matrix C
		C = (double **) malloc(A_rows * sizeof(double *));
		for(i=0; i<A_rows ;i++){
			C[i] = (double *) malloc(B_columns * sizeof(double));
		}
	} 

	// send a block to each process
	if(rank == 0) {
		int i;
		for(i = 1; i < size; i++){
			// Each process gets one A and one B block
			MPI_Send((A_array + (i * A_local_block_size)), A_local_block_size, MPI_DOUBLE, i, 0, cartesian_grid_communicator);
			MPI_Send((B_array + (i * B_local_block_size)), B_local_block_size, MPI_DOUBLE, i, 0, cartesian_grid_communicator);
		}
		// Filling in block 0
		for(i = 0; i < A_local_block_size; i++){
			A_local_block[i] = A_array[i];
		}
		for(i = 0; i < B_local_block_size; i++){
			B_local_block[i] = B_array[i];
		}
	} else {
		MPI_Recv(A_local_block, A_local_block_size, MPI_DOUBLE, 0, 0, cartesian_grid_communicator, &status);
		MPI_Recv(B_local_block, B_local_block_size, MPI_DOUBLE, 0, 0, cartesian_grid_communicator, &status);
	}

	// cannon's algorithm
	double compute_time = 0, mpi_time = 0, start;
	//start = MPI_Wtime();
	
	MPI_Info info;
	int info_res = MPI_Info_create(&info);
	//int MPI_Info_set(MPI_Info info, const char *key, const char *value)
	char* no_locks = "no_locks";
	char* true_msg = "true";
	char* same_size = "same_size";	
	char* same_disp_unit = "same_disp_unit";	
	MPI_Info_set(info, no_locks, true_msg);
	MPI_Info_set(info, same_size, true_msg);
	MPI_Info_set(info, same_disp_unit, true_msg);

	MPI_Group row_group;
	MPI_Group col_group;
	MPI_Comm_group(row_communicator, &row_group);
	MPI_Comm_group(column_communicator, &col_group);
	
	MPI_Win_create(a_block, 2 * A_local_block_size, sizeof(double), MPI_INFO_NULL, row_communicator, &win_A);
	MPI_Win_create(b_block, 2 * B_local_block_size, sizeof(double), MPI_INFO_NULL, column_communicator, &win_B);
	//MPI_Win_create(A_local_block_new, A_local_block_size, sizeof(double), MPI_INFO_NULL, row_communicator, &win_A);
	//MPI_Win_create(B_local_block_new, B_local_block_size, sizeof(double), MPI_INFO_NULL, column_communicator, &win_B);
	
	MPI_Win_set_info(win_A, info);
	MPI_Win_set_info(win_B, info);
	
	//MPI_Win_fence(0, win_A);
	//MPI_Win_fence(0, win_B);
	//mpi_time += MPI_Wtime() - start;
	int cannon_block_cycle;
	double put_time = 0;
	double fence_time = 0;
	double wt;
	int C_index, A_row, A_column, B_column;
	int displacement_A = A_local_block_size;
	int displacement_B = B_local_block_size;
	int current_disp_A = 0;
	int current_disp_B = 0;
	for(cannon_block_cycle = 0; cannon_block_cycle < sqrt_size; cannon_block_cycle++){
		// compute partial result for this block cycle
		start = MPI_Wtime();
		
		MPI_Win_post(row_group, 0, win_A);
		MPI_Win_start(row_group, 0, win_A);
		MPI_Put(A_local_block, A_local_block_size, MPI_DOUBLE, (coordinates[1] + sqrt_size - 1) % sqrt_size, current_disp_A, A_local_block_size, MPI_DOUBLE, win_A);
		
		
		MPI_Win_post(col_group, 0, win_B);
		MPI_Win_start(col_group, 0, win_B);
		MPI_Put(B_local_block, B_local_block_size, MPI_DOUBLE, (coordinates[0] + sqrt_size - 1) % sqrt_size, current_disp_B, B_local_block_size, MPI_DOUBLE, win_B);
		
		
		
		wt = MPI_Wtime() - start;
		put_time += wt;
		mpi_time += wt;
		
		start = MPI_Wtime();
		for(C_index = 0, A_row = 0; A_row < A_local_block_rows; A_row++){
			for(B_column = 0; B_column < B_local_block_columns; B_column++, C_index++){
				for(A_column = 0; A_column < A_local_block_columns; A_column++){
					C_local_block[C_index] += A_local_block[A_row * A_local_block_columns + A_column] *
						B_local_block[A_column * B_local_block_columns + B_column];
				}
			}
		}
		
	
		compute_time += MPI_Wtime() - start;
		start = MPI_Wtime();
		
		MPI_Win_complete(win_A);
		MPI_Win_complete(win_B);
		MPI_Win_wait(win_A);
		MPI_Win_wait(win_B);
				
		wt = MPI_Wtime() - start;
		fence_time += wt;
		mpi_time += wt;
		double * tmp = A_local_block;
		A_local_block = A_local_block_new;
		A_local_block_new = tmp;

		tmp = B_local_block;
		B_local_block = B_local_block_new;
		B_local_block_new = tmp;
		
		current_disp_A ^= displacement_A;
		current_disp_B ^= displacement_B;
	}
 	//MPI_Win_Free(win_A);
	//MPI_Win_Free(win_B);
	double overall_mpi_time;
	double overall_compute_time;
	double opt, oft;
	MPI_Reduce(&mpi_time, &overall_mpi_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&put_time, &opt, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&fence_time, &oft, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&compute_time, &overall_compute_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	// get C parts from other processes at rank 0
	if(rank == 0) {
		for(i = 0; i < A_local_block_rows * B_local_block_columns; i++){
			C_array[i] = C_local_block[i];
		}
		int i;
		for(i = 1; i < size; i++){
			MPI_Recv(C_array + (i * A_local_block_rows * B_local_block_columns), A_local_block_rows * B_local_block_columns, 
				MPI_DOUBLE, i, 0, cartesian_grid_communicator, &status);
		}
	} else {
		MPI_Send(C_local_block, A_local_block_rows * B_local_block_columns, MPI_DOUBLE, 0, 0, cartesian_grid_communicator);
	}

	// generating output at rank 0
	if (rank == 0) {
		// convert the ID array into the actual C matrix 
		int i, j, k, row, column;
		for (i = 0; i < sqrt_size; i++){  // block row index
			for (j = 0; j < sqrt_size; j++){ // block column index
				for (row = 0; row < A_local_block_rows; row++){
					for (column = 0; column < B_local_block_columns; column++){
						C[i * A_local_block_rows + row] [j * B_local_block_columns + column] = 
							C_array[((i * sqrt_size + j) * A_local_block_rows * B_local_block_columns) 
							+ (row * B_local_block_columns) + column];
					}
				}
			}
		}

		printf("(%d,%d)x(%d,%d)=(%d,%d)\n", A_rows, A_columns, B_rows, B_columns, A_rows, B_columns);
		printf("Computation time: %lf\n", overall_compute_time/size);
		printf("MPI time:         %lf\n", overall_mpi_time/size);
		printf("Fence time:         %lf\n", oft/size);
		printf("Put time:         %lf\n", opt/size);				

		if (argc == 4){
			// present results on the screen
			/*printf("\nA( %d x %d ):\n", A_rows, A_columns);
			for(row = 0; row < A_rows; row++) {
				for(column = 0; column < A_columns; column++)
					printf ("%7.3f ", A[row][column]);
				printf ("\n");
			}
			printf("\nB( %d x %d ):\n", B_rows, B_columns);
			for(row = 0; row < B_rows; row++){
				for(column = 0; column < B_columns; column++)
					printf("%7.3f ", B[row][column]);
				printf("\n");
			}
			printf("\nC( %d x %d ) = AxB:\n", A_rows, B_columns);
			for(row = 0; row < A_rows; row++){
				for(column = 0; column < B_columns; column++)
					printf("%7.3f ",C[row][column]);
				printf("\n");
			}
*/

			printf("\nPerforming serial consistency check. Be patient...\n");
			fflush(stdout);
			int pass = 1;
			double temp;
			for(i=0; i<A_rows; i++){
				for(j=0; j<B_columns; j++){
					temp = 0;
					for(k=0; k<B_rows; k++){
						temp += A[i][k] * B[k][j];
					}
					//printf("%7.3f ", temp);
					if(temp != C[i][j]){
						pass = 0;
					}
				}
				//printf("\n");
			}
			if (pass) printf("Consistency check: PASS\n");
			else printf("Consistency check: FAIL\n");
		}	
	}

	// free all memory
	if(rank == 0){
		int i;
		for(i = 0; i < A_rows; i++){
			free(A[i]);
		}
		for(i = 0; i < B_rows; i++){
			free(B[i]);
		}
		for(i = 0; i < A_rows; i++){
			free(C[i]);
		}
		free(A);
		free(B);
		free(C);
		free(A_array);
		free(B_array);
		free(C_array);
	}
	free(A_local_block);
	free(B_local_block);
	free(C_local_block);

	// finalize MPI
	MPI_Finalize();
}


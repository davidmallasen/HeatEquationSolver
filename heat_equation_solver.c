#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define DEBUG 1

// Length of the 2D domain
#define L 4//3360 // Multiple of 4, 6, 8, 10, 12, 14, 16 (up to 256 MPI processes)

// Width of the grid spacing
#define H 1

// Thermal diffusivity of the medium
#define K 2

// Fixed boundary temperature
#define TEMP_BOUND 100

// Starting center temperature
#define TEMP_CENTER 0

// End time of of the simulation (starting at 0 and with TIME_DELTA increments)
#define TIME_END 0.2

// Delta of time between iterations
#define TIME_DELTA 0.1


void init_temp(double *mat, int dim, int west_neigh, int east_neigh,
               int north_neigh, int south_neigh);
void update_steps(double *prev_mat, double *mat, int dim);
void halo_exchange(double *mat, int dim, int west_neigh, int east_neigh, 
                   int north_neigh, int south_neigh, MPI_Comm cart_comm, 
                   MPI_Datatype mpi_column_t);
void print_matrix(double *local_mat, int local_dim, int mat_dim, 
                  int grid_side_len, int rank, MPI_Comm cart_comm);


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int grid_side_len = (int) sqrt((double) size);

    int local_dim = L / grid_side_len;
    if (L % grid_side_len != 0) {
        printf("Error. L must be a multiple of grid_side_len\n");
        return 1;
    }

    // Create cartesian topology
    int dim_sizes[2];
    int wrap_around[2];
    MPI_Comm cart_comm;
    int rank;

    dim_sizes[0] = grid_side_len;
    dim_sizes[1] = grid_side_len;
    wrap_around[0] = 0;
    wrap_around[1] = 0;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes, wrap_around, 1, &cart_comm);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Identify neighbors
    int west_neigh, east_neigh, north_neigh, south_neigh;

    MPI_Cart_shift(cart_comm, 0, 1, &west_neigh, &east_neigh);
    MPI_Cart_shift(cart_comm, 1, 1, &south_neigh, &north_neigh);

    // Declare column type
    MPI_Datatype mpi_column_t;

    MPI_Type_vector(local_dim, 1, local_dim + 2, MPI_DOUBLE, &mpi_column_t);
    MPI_Type_commit(&mpi_column_t);

    // Allocate and initialize the matrix (including ghost cells)
    double *mat = malloc((local_dim + 2) * (local_dim + 2) * sizeof(double));
    double *next_mat = malloc((local_dim + 2) * (local_dim + 2) * sizeof(double));

    init_temp(mat, local_dim, west_neigh, east_neigh, north_neigh, south_neigh);

    // Computation
    double start_time, stop_time, elapsed_time;
    start_time = MPI_Wtime();

    for (float t = 0; t < TIME_END; t += TIME_DELTA) {
        halo_exchange(mat, local_dim, west_neigh, east_neigh, north_neigh, 
                      south_neigh, cart_comm, mpi_column_t);
        update_steps(mat, next_mat, local_dim);
    }

    stop_time = MPI_Wtime();
    elapsed_time = stop_time - start_time;

    // Print results
    #if DEBUG
        print_matrix(mat, local_dim, L, grid_side_len, rank, cart_comm);
    #else
        if (rank == 0) {
            printf("Total time: %f\n", elapsed_time);
        }
    #endif

    // Free resources and finalize
    free(mat);
    free(next_mat);

    MPI_Type_free(&mpi_column_t);

    MPI_Finalize();
    
    return 0;
}

void init_temp(double *mat, int dim, int west_neigh, int east_neigh, 
               int north_neigh, int south_neigh) {
    double west_value = west_neigh == MPI_PROC_NULL ? TEMP_BOUND : TEMP_CENTER;
    double east_value = east_neigh == MPI_PROC_NULL ? TEMP_BOUND : TEMP_CENTER;
    double north_value = north_neigh == MPI_PROC_NULL ? TEMP_BOUND : TEMP_CENTER;
    double south_value = south_neigh == MPI_PROC_NULL ? TEMP_BOUND : TEMP_CENTER;

    for (int i = 0; i < dim + 2; i++) {
        mat[i * (dim + 2)] = west_value;
        mat[i * (dim + 2) + dim + 1] = east_value;
        mat[i] = north_value;
        mat[(dim + 1) * (dim + 2) + i] = south_value;
    }

    for (int i = 1; i < dim + 1; i++) {
        for (int j = 1; j < dim + 1; j++) {
            mat[i * (dim + 2) + j] = TEMP_CENTER;
        }
    }
}

void update_steps(double *mat, double *next_mat, int dim) {
    double tmp = (K * TIME_DELTA);
    double h2 = H * H;
    int dim2 = dim + 2;
    for (int i = 1; i < dim + 1; i++) {
        for (int j = 1; j < dim + 1; j++) {
            next_mat[i * dim2 + j] = mat[i * dim2 + j] + tmp * (
                (mat[(i + 1) * dim2 + j] - 2 * mat[i * dim2 + j] + mat[(i - 1) * dim2 + j]) / h2 + 
                (mat[i * dim2 + (j + 1)] - 2 * mat[i * dim2 + j] + mat[i * dim2 + (j - 1)]) / h2
            );
        }
    }

    for (int i = 1; i < dim + 1; i++) { // Borders will be updated from neighbors
        for (int j = 1; j < dim + 1; j++) {
            mat[i * dim2 + j] = next_mat[i * dim2 + j];
        }
    }
}

void halo_exchange(double *mat, int dim, int west_neigh, int east_neigh, 
                   int north_neigh, int south_neigh, MPI_Comm cart_comm, 
                   MPI_Datatype mpi_column_t) {
    // Send west border (second column) and receive east border (last column)
    MPI_Sendrecv(&mat[(dim + 2) + 1], 1, mpi_column_t, west_neigh, 0, 
                 &mat[(dim + 2) + dim + 1], 1, mpi_column_t, east_neigh, 0, 
                 cart_comm, MPI_STATUS_IGNORE);

    // Send east border (second to last column) and receive west border (first column)
    MPI_Sendrecv(&mat[(dim + 2) + dim], 1, mpi_column_t, east_neigh, 1, 
                 &mat[(dim + 2)], 1, mpi_column_t, west_neigh, 1, 
                 cart_comm, MPI_STATUS_IGNORE);

    // Send north border (second row) and receive south border (last row)
    MPI_Sendrecv(&mat[(dim + 2) + 1], dim, MPI_DOUBLE, north_neigh, 2, 
                 &mat[(dim + 1) * (dim + 2) + 1], dim, MPI_DOUBLE, south_neigh, 2, 
                 cart_comm, MPI_STATUS_IGNORE);

    // Send south border (second to last row) and receive north border (first row)
    MPI_Sendrecv(&mat[dim * (dim + 2) + 1], dim, MPI_DOUBLE, south_neigh, 3, 
                 &mat[1], dim, MPI_DOUBLE, north_neigh, 3, 
                 cart_comm, MPI_STATUS_IGNORE);
}

void print_matrix(double *local_mat, int local_dim, int mat_dim, int grid_side_len, int rank, MPI_Comm cart_comm) {
    int source;
    int coords[2];

    // Process 0 receives all the matrix blocks and writes them
    if (rank == 0) {
        double* temp_row = (double*) malloc(local_dim * sizeof(double));

        for (int mat_row = 0; mat_row < mat_dim; mat_row++) {
            // Local row corresponding to current column pointer
            coords[0] = mat_row / local_dim;

            for (int grid_col = 0; grid_col < mat_dim; grid_col += local_dim) {
                coords[1] = grid_col / local_dim;
                source = coords[1] * grid_side_len + grid_side_len - coords[0] - 1;
                if (source == 0) {
                    for(int mat_col = 0; mat_col < local_dim; mat_col++) {
                        printf("%04.2f ", local_mat[(mat_row % local_dim + 1) * (local_dim + 2) + (mat_col + 1)]);
                    }
                }
                else {
                    MPI_Recv(temp_row, local_dim, MPI_DOUBLE, source, 0, cart_comm,
                             MPI_STATUS_IGNORE);
                    for(int mat_col = 0; mat_col < local_dim; mat_col++) {
                        printf("%04.2lf ", temp_row[mat_col]);
                    }
                }
            }
            printf("\n");
        }

        free(temp_row);
    } 
    else { // Send matrix block row by row to process 0
        for (int mat_row = 1; mat_row < local_dim + 1; mat_row++) {
            MPI_Send(&local_mat[mat_row * (local_dim + 2) + 1], local_dim, MPI_DOUBLE, 0, 0, cart_comm);
        }
    }
}

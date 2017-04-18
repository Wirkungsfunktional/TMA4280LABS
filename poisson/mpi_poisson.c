#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>


#define PI 3.14159265358979323846
#define true 1
#define false 0
#define CHECK(x) if (x==NULL) {printf("Malloc error!!!");exit(0);}
#define IO_OPT 1

typedef double real;
typedef int bool;

void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);



real rhs(real x, real y) {
    return exp(x)*sin(2*PI*x)*sin(2*PI*y);
    //return 2 * (y - y*y + x - x*x);
}

real* make_1D(int n) {
    return (real*) malloc(n*sizeof(real));
    }
    
real** make_2D(int n, int m) {
    real** a = (real**) malloc(n * sizeof(real*));
    a[0] = (real*) malloc(n * m * sizeof(real)); 
    for (int i=1;i<n;i++) a[i] = a[i-1] + m;
    return a;
    }

int main(int argc, char **argv)
{
    int size, rank, tag, i;
    double ta, te;  
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ta = MPI_Wtime(); 
    
    int n = atoi(argv[1]);
    int m = n - 1;
    real h = 1.0 / n;
    int num_row = n / size;
    int last_num_row = m - num_row * (size-1);
    int local_num_row;
    
    if (last_num_row > num_row) {
        num_row += 1;
        last_num_row = m - num_row * (size-1);
    }
    
    if (rank == size - 1) {
        local_num_row = last_num_row;
    } else {
        local_num_row = num_row;
    }
    
    if (rank==0) {
        printf("%d\n", local_num_row);
        printf("%d\n", last_num_row);
    }


    int scounts[size];CHECK(scounts);
    int sdispls[size];CHECK(sdispls);
    int rcounts[size];CHECK(rcounts);
    int rdispls[size];CHECK(rdispls);
    
    
    real *b_row = make_1D(m); CHECK(b_row);
    real *sendbuffer = make_1D(size*num_row*num_row); CHECK(sendbuffer);
    real *recvbuffer = make_1D(size*num_row*num_row); CHECK(recvbuffer);
    int nn = 4*n;
    real *z = make_1D(nn); CHECK(z);
    for (int i=0 ; i< size; i++) {
        scounts[i] = num_row*num_row;
        sdispls[i] = i*num_row*num_row;
        rcounts[i] = num_row*num_row;
        rdispls[i] = i*num_row*num_row;
    }

    real *diag = make_1D(m); CHECK(diag);
    real *grid = make_1D(n+1); CHECK(grid);
    #pragma omp parallel for
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
        grid[i] = i * h;
    }
    grid[m] = m * h;
    grid[n] = n * h;
    
    
    for (size_t i = 0; i < local_num_row; i++) {
        #pragma omp parallel for
        for (size_t j = 0; j < m; j++) {
            b_row[j] = h * h * rhs(grid[(num_row*rank + i)+1], grid[j+1]);
        }
        fst_(b_row, &n, z, &nn);
        for (int k=0;k<size-1;k++) {
            for (int j=0; j<num_row; j++) {
                sendbuffer[k*num_row*num_row + i*num_row + j] = b_row[k*num_row + j];
            }
        }
        #pragma omp parallel for
        for (int j=0; j<last_num_row; j++) {
            sendbuffer[(size-1)*num_row*num_row + i*num_row + j] = b_row[(size-1)*num_row + j];
        }
    }
    MPI_Alltoallv(sendbuffer, scounts, sdispls, MPI_DOUBLE, recvbuffer, rcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
    
    for (int i=0; i < local_num_row;i++) {
        #pragma omp parallel for
        for (int j=0; j < m;j++) {
            b_row[j] = recvbuffer[j*num_row + i];
        }
        fstinv_(b_row, &n, z, &nn);
        #pragma omp parallel for
        for (size_t j = 0; j < m; j++) {
            b_row[j] = b_row[j] / (diag[(num_row*rank + i)] + diag[j]);
        }
        fst_(b_row, &n, z, &nn);
        for (int k=0;k<size-1;k++) {
            for (int j=0; j<num_row; j++) {
                sendbuffer[k*num_row*num_row + i*num_row + j] = b_row[k*num_row + j];
            }
        }
        #pragma omp parallel for
        for (int j=0; j<last_num_row; j++) {
            sendbuffer[(size-1)*num_row*num_row + i*num_row + j] = b_row[(size-1)*num_row + j];
        }
    }
    MPI_Alltoallv(sendbuffer, scounts, sdispls, MPI_DOUBLE, recvbuffer, rcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
    free(sendbuffer);
    free(b_row);
    
    real **b = make_2D(local_num_row, m);CHECK(b);
    for (int i=0; i < local_num_row;i++) {
        #pragma omp parallel for
        for (int j=0; j < m;j++) {
            b[i][j] = recvbuffer[j*num_row + i];
        }
        fstinv_(b[i], &n, z, &nn);
    }
    free(z);
    free(recvbuffer);
    free(grid);
    free(diag);
    
    if (IO_OPT) {
        MPI_File fh;
        MPI_Offset my_offset = sizeof(real) * rank * m * (num_row);
        int number_of_data = m * local_num_row;
        MPI_File_open(MPI_COMM_WORLD, "test.txt", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
        MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
        for (int i=0; i< local_num_row; i++)
            MPI_File_write(fh, b[i], m, MPI_DOUBLE, NULL);
        MPI_File_close(&fh);
    }
    
    te = MPI_Wtime(); 
    
    real loc_max, global_max;
    double t_max, t;
    
    t = te -ta;
    MPI_Allreduce(&t, &t_max, 1, MPI_DOUBLE, MPI_MAX,
                            MPI_COMM_WORLD);

            
    
    loc_max = 0.0;
    for (size_t i = 0; i < local_num_row; i++) {
        for (size_t j = 0; j < m; j++) {
            loc_max = loc_max > b[i][j] ? loc_max : b[i][j];
        }
    }
    MPI_Allreduce(&loc_max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
                        MPI_COMM_WORLD);
    //free(b); // WARN
    if (rank == 0) {
        printf("%e\n", global_max - 1/(5*PI*PI));
        printf("%0.7f\n", t_max);
    }

    
    MPI_Finalize();
    
    
	
	return 0;
}


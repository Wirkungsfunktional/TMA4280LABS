#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define PI 3.14159265358979323846
#define true 1
#define false 0

typedef double real;
typedef int bool;

void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

real rhs(real x, real y) {
    return sin(PI*x)*sin(2*PI*y);//2 * (y - y*y + x - x*x);
}

int main(int argc, char **argv)
{
    int size, rank, tag, i;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int n = atoi(argv[1]);
    int m = n - 1;
    real h = 1.0 / n;
    int num_row = n / size;
    int last_num_row = m - num_row * (size-1);
    int local_num_row;
    
    if (rank == size - 1) {
        local_num_row = last_num_row;
    } else {
        local_num_row = num_row;
    }
    //printf("%d", last_num_row);
    
    
    int scounts[size*sizeof(int)];
    int sdispls[size*sizeof(int)];
    int rcounts[size*sizeof(int)];
    int rdispls[size*sizeof(int)];
    

    
    //printf("%i\n", local_num_row);
    real grid[n+1];
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }
    
    real diag[m];
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    
    real b[local_num_row][m];
    real b_row[m];
    real sendbuffer[size][local_num_row][num_row];
    real recvbuffer[size*num_row*num_row];
    int nn = 4*n;
    real z[nn];
    
    
    //printf("rtn1\n");
    for (size_t i = 0; i < local_num_row; i++) {
        for (size_t j = 0; j < m; j++) {
            b_row[j] = h * h * rhs(grid[(num_row*rank + i)+1], grid[j+1]);
        }
        //fst_(b[i], &n, z, &nn);
        for (int k=0;k<size;k++) {
            for (int j=0; j<num_row; j++) {
                sendbuffer[k][i][j] = b_row[k*num_row + j];
            }
        }
        for (int j=0; j<num_row; j++) {
            sendbuffer[(size-1)][i][j] = b_row[(size-1)*num_row + j];
        }
    }
    for (int i=0 ; i< size; i++) {
        scounts[i] = local_num_row*num_row;
        sdispls[i] = i*local_num_row*num_row;
        rcounts[i] = local_num_row*num_row;
        rdispls[i] = i*num_row*num_row;
    }
    scounts[size-1] = last_num_row*num_row;
    //rcounts[size-1] = last_num_row*num_row;
    
    MPI_Alltoallv(sendbuffer, scounts, sdispls, MPI_DOUBLE, recvbuffer, rcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
    for (int i=0; i<num_row; i++) {
        for (int j=0;j<m;j++) {
            b[i][j] = recvbuffer[i + num_row*j];
        }
    }
    
    
    
    
    
    /*
    printf("rtn2\n");
    for (size_t i = 0; i < num_row; i++) {
        fstinv_(bt[i], &n, z, &nn);
        for (size_t j = 0; j < m; j++) {
            bt[i][j] = bt[i][j] / (diag[(num_row*rank + i)] + diag[j]);
        }
        fst_(bt[i], &n, z, &nn);
        MPI_Alltoallv(bt[i], scounts, sdispls, MPI_DOUBLE, storage, rcounts, rdispls, MPI_DOUBLE, MPI_COMM_WORLD);
        for (int k=0;k<size;k++) {
            for (int l=0;l<num_row;l++) {
                b[l][k*num_row + i] = storage[k*num_row + l];
            }
        }
    }
    for (size_t i = 0; i < num_row; i++) {
        fstinv_(b[i], &n, z, &nn);
    }
    printf("rtn3\n");
    real loc_max, global_max;
    loc_max = 0.0;
    for (size_t i = 0; i < num_row; i++) {
        for (size_t j = 0; j < m; j++) {
            loc_max = loc_max > b[i][j] ? loc_max : b[i][j];
        }
    }
    
    
    MPI_Allreduce(&loc_max, &global_max, 1, MPI_DOUBLE, MPI_MAX,
              MPI_COMM_WORLD);
    
    
    */
    num_row = n / size;
    
    MPI_File fh;
    MPI_Offset my_offset = sizeof(real) * rank * m * (num_row);
    int number_of_data = m * local_num_row;
    MPI_File_open(MPI_COMM_WORLD, "test.txt", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
    
    MPI_File_seek(fh, my_offset, MPI_SEEK_SET);
    MPI_File_write(fh, b, number_of_data, MPI_DOUBLE, NULL);
    MPI_File_close(&fh);
    
        
    if (rank == 0) {
        //printf("%f", global_max);
    }
        

    
    
    
    MPI_Finalize();
    
    
	
	return 0;
}


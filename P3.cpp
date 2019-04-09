#include "mmreader.hpp"
#include <time.h>
#include <iostream>
#include <sys/time.h>
#include <unistd.h>

#define NUM_THREADS 12

bool
SCsrMatrixfromFile(struct sparse_mtx *A, const char* filePath)
{
    // Check that the file format is matrix market; the only format we can read right now
    // This is not a complete solution, and fails for directories with file names etc...
    // TODO: Should we use boost filesystem?
    std::string strPath( filePath );
    if( strPath.find_last_of( '.' ) != std::string::npos )
    {
        std::string ext = strPath.substr( strPath.find_last_of( '.' ) + 1 );
        if( ext != "mtx" )
        {
            std::cout << "Reading file name error" << std::endl;
            return false;
        }
    }
    else
        return false;

    // Read data from a file on disk into buffers
    // Data is read natively as COO format with the reader
    MatrixMarketReader mm_reader;
    if( mm_reader.MMReadFormat(filePath) )
        return false;

    // JPA: Shouldn't that just be an assertion check? It seems to me that
    // the user have to call clsparseHeaderfromFile before calling this function,
    // otherwise the whole pCsrMatrix will be broken;
    A->nrow = mm_reader.GetNumRows( );
    A->ncol = mm_reader.GetNumCols( );
    A->nnze = mm_reader.GetNumNonZeroes( );

    A->row = (int32_t *)malloc((A->nrow + 1) * sizeof(int32_t));
    A->val = (float *)malloc(A->nnze * sizeof(float));
    A->col = (int32_t *)malloc(A->nnze * sizeof(int32_t));

    if(A->row == NULL || A->col == NULL || A->val == NULL)
    {
        if(A->row == NULL)
            free((void *)A->row);
        if(A->col == NULL)
            free((void *)A->col);
        if(A->val == NULL)
            free((void *)A->val);
        return false;
    }

    //  The following section of code converts the sparse format from COO to CSR
    Coordinate* coords = mm_reader.GetUnsymCoordinates( );

    std::sort( coords, coords + A->nnze, CoordinateCompare );

    int32_t current_row = 1;

    A->row[ 0 ] = 0;

    for (int32_t i = 0; i < (int32_t)A->nnze; i++)
    {
        A->col[ i ] = coords[ i ].y;
        A->val[ i ] = coords[ i ].val;

        while( coords[ i ].x >= current_row )
            A->row[ current_row++ ] = i;
    }

    A->row[ current_row ] = A->nnze;

    while( current_row <= (int32_t)A->nrow )
        A->row[ current_row++ ] = A->nnze;

    return true;
}

void multiply_single(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C)
{
    // TODO: Implement matrix multiplication with single thread. C=A*B
    uint32_t col_idx = 0;
    uint32_t i, j;
    for (i=0; i<A->nrow-1; i++){
        for (j=A->row[i];j<(uint32_t)A->row[i+1];j++){
            col_idx = A->col[j];
            C->val[i*(C->ncol)+col_idx] += (A->val[col_idx])*(B->val[col_idx*(B->ncol)+i]);
        }
    }
    // last row handling
    for (j=A->row[i];j<(uint32_t)A->ncol;j++){
        col_idx = A->col[j];
        C->val[i*(C->ncol)+col_idx] += (A->val[col_idx])*(B->val[col_idx*(B->ncol)+i]);
    }
}

void multiply_pthread(struct sparse_mtx *A, struct dense_mtx *B, struct dense_mtx *C)
{
    // TODO: Implement matrix multiplication with pthread. C=A*B
    uint32_t col_idx = 0;
    uint32_t i, j;
     
}

void init_zeros(struct dense_mtx *A, uint32_t nrow, uint32_t ncol)
{
    A->nrow = nrow;
    A->ncol = ncol;
    A->val = (float *)malloc(sizeof(float)*nrow*ncol);
    for (uint32_t i=0; i<nrow*ncol; i++){
        A->val[i] = 0;
    }
}

double time_elapse(struct timespec start, struct timespec end) {
    double time = end.tv_sec - start.tv_sec;
    time = time*1000000 + (double)(end.tv_nsec - start.tv_nsec)/1000;
    return time;
}

int main(int argc, char **argv)
{
    struct sparse_mtx A;
    if(!SCsrMatrixfromFile(&A, argv[1]))
    {
        std::cout << "read failed." << std::endl;
        return 0;
    }

    struct dense_mtx B;
    B.nrow = A.ncol;
    B.ncol = atoi(argv[2]);
    if(B.ncol < 0)
    {
        free(A.row);
        free(A.col);
        free(A.val);
        std::cerr << "Invalid argument for the number of columns of B." << std::endl;
        return 0;
    }
    B.val = (float *)malloc(sizeof(float) * B.nrow * B.ncol);

    srand((unsigned int)time(NULL));
    for(uint32_t i = 0; i < B.nrow; i++)
    {
        for(uint32_t j = 0; j < B.ncol; j++)
        {
            B.val[B.ncol * i + j] = ((float)rand()/(float)(RAND_MAX)) * ((rand() % 2) ? 1.0f : -1.0f);
        }
    }

    struct dense_mtx C1, C2;
    init_zeros(&C1, A.nrow, B.ncol);
    init_zeros(&C2, A.nrow, B.ncol);

    std::cout << "Single Thread Computation Start" << std::endl;
    
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    multiply_single(&A, &B, &C1);
    clock_gettime(CLOCK_MONOTONIC, &end);
    std::cout << "Single Thread Computation End: " << time_elapse(start, end)  << " us." << std::endl;
    std::cout << "Multi Thread Computation Start" << std::endl;
    clock_gettime(CLOCK_MONOTONIC, &start);
    multiply_pthread(&A, &B, &C2);
    clock_gettime(CLOCK_MONOTONIC, &end);
    std::cout << "Multi Thread Computation End: " << time_elapse(start, end) << " us." << std::endl;

    // TODO: Testing Code by comparing C1 and C2

    free(A.row);
    free(A.col);
    free(A.val);
    free(B.val);
    if(C1.val != NULL)
        free(C1.val);
    if(C2.val != NULL)
        free(C2.val);
    
    return 0;
}

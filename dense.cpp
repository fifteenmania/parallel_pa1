#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#define NUM_THREADS 6
#define MEQ_THRES 0.0001;

using namespace std;

typedef struct{
    double *A;
    double *B;
    double *C;
    int n;
    int pid;
}mm_args;

double *init_zeros(int n)
{
    int ndim = n*n;
    double *C = (double*) malloc(ndim*sizeof(double));
    if (!C){
        return NULL;
    }
    for (int i=0; i<ndim; i++){
        C[i] = 0;
    }
    return C;
}

double *init_rand(int n)
{
    int ndim = n*n;
    double *A = (double*) malloc(ndim*sizeof(double));
    if (!A){
        return NULL;
    }
    srand48(10);
    for (int i=0; i<ndim; i++){
        A[i] = drand48();
    }
    return A;
}

void mmult_naive(double *A, double *B, double *C, int n)
{
    int ndim = n*n;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                C[i*n+j] += A[i*n+k]*B[n*k+j];
            }
        }
    }
    return;
}

void _mmult_naive_par(void *args_pasd)
{
    mm_args* args = (mm_args *) args_pasd;
    double *A = (double *)args->A;
    double *B = (double *)args->B;
    double *C = (double *)args->C;
    int pid = (int) args->pid;
    int n = (int) args->n;
    int st_part = (int)((double)pid*(double)n/NUM_THREADS);
    int ed_part = (int)((double)(pid+1)*(double)n/NUM_THREADS);
    for (int i=st_part; i<ed_part; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                C[i*n+j] += A[i*n+k]*B[n*k+j];
            }
        }
    }
    return;
}

void mmult_naive_par(double *A, double *B, double *C, int n)
{
    pthread_t p_threads[NUM_THREADS];
    mm_args args[NUM_THREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    for (int i=0; i<NUM_THREADS; i++){
        args[i].A = A;
        args[i].B = B;
        args[i].C = C;
        args[i].n = n;
        args[i].pid = i;
        pthread_create(&p_threads[i], &attr, _mmult_naive_par, (void *) &(args[i]));
    }
    return;
}

double mres(double *A, double *B, int n)
{
    double l2n=0;
    double res;
    int ndim = n*n;
    for (int i=0; i<ndim; i++){
        res = A[i] - B[i];
        l2n = res*res;
    }
    return l2n;
}

bool mequal(double *A, double *B, int n)
{
    return mres(A, B, n) < MEQ_THRES;
}

void _print_mat(double *A, int n)
{
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            cout << A[n*i+j] << " ";
        }
        cout << endl;
    }
    return;
}

void _print_mats(double *A, double *B, double *C, int n)
{
    cout << "Matrix A: " << endl;
    _print_mat(A, n);
    cout << "Matrix B: " << endl;
    _print_mat(B, n);
    cout << "Matrix C: " << endl;
    _print_mat(C, n);
    return;
}


int main(int argc, char **argv)
{
    int n = 1000;
    double *A, *B;
    double *C1, *C2;
    static struct timespec begin1, end1, begin2, end2;
    double bench_time;
    // initialization
    A = init_rand(n);
    B = init_rand(n);
    C1 = init_zeros(n);
    C2 = init_zeros(n);
    // multiplication routine
    clock_gettime(CLOCK_MONOTONIC, &begin1);
    mmult_naive_par(A, B, C1, n);
    clock_gettime(CLOCK_MONOTONIC, &end1);
    // print benchmark results
    bench_time = (double) (end1.tv_nsec - begin1.tv_nsec);
    bench_time = bench_time/100000000000 + end1.tv_sec - begin1.tv_sec;
    cout << "Time elapsed: " << bench_time << " sec" << endl;
    free(A);
    free(B);
    free(C);
    return 0;
}


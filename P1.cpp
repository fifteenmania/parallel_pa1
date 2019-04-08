#include "Matrix.h"
#include <iostream>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

double time_elapse(struct timespec begin, struct timespec end)
{
    double time = end.tv_sec - begin.tv_sec;
    time = time*1000 + (double)(end.tv_nsec - begin.tv_nsec)/1000000;
    return time;
}

int main(int argc, char **argv)
{
    int n, p;
    Matrix A, B, C;
    Matrix D;
    if (argc !=3){
        cout << "Usage: ./P1 [n]:size of matrix [p]:number of threads" << endl;
        return 0;
    }
    n = atoi(argv[1]);
    p = atoi(argv[2]);
    if (n<1 || n>100000){
        cout << "Invalid input" << endl;
        return 0;
    }
    if (p<1 || p>1024){
        cout << "Invalid input" << endl;
        return 0;
    }
    A.init_rand(n, 10);
    B.init_rand(n, 20);
    C.init_zeros(n);
    D.init_zeros(n);

    struct timespec begin, end;
    //clock_gettime(CLOCK_MONOTONIC, &begin);
    C.set_mmult_ordered(A, B);
    //clock_gettime(CLOCK_MONOTONIC, &end);
    //cout << "Single:   " << time_elapse(begin, end) << " ms" << endl;
    
    clock_gettime(CLOCK_MONOTONIC, &begin);
    D.set_mmult_par(A, B, p);
    clock_gettime(CLOCK_MONOTONIC, &end);
    cout << "Parallel: " << time_elapse(begin, end) << " ms" << endl;
    cout << "Correctness: " << C.is_equal(D) << endl;
    return 0;
}

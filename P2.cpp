#include "Matrix.h"
#include <iostream>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>

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
    Matrix A, A0;
    Vector b, b0;
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
    if (p<1 || p>512){
        cout << "Invalid input" << endl;
        return 0;
    }
    A.init_rand(n, 0);//time(NULL));
    A0 = A;
    b.init_rand(n, 1);//time(NULL));
    b0 = b;

    struct timespec begin, end;
    //clock_gettime(CLOCK_MONOTONIC, &begin);
    //clock_gettime(CLOCK_MONOTONIC, &end);
    //cout << "Single:      " << time_elapse(begin, end) << " ms" << endl;
    //A.set_gauss_elim_par(&b, p);
    clock_gettime(CLOCK_MONOTONIC, &begin);
    A.set_gauss_elim_par(&b, p);
    A.set_backsub(&b, p);
    clock_gettime(CLOCK_MONOTONIC, &end);
    cout << "Parallel:    " << time_elapse(begin, end) << " ms" << endl;
    
    b = b.vmult(&A0);
    double resid = b.residual_diff(b0);
    cout << "Residual:    " << resid << endl;
    cout << "Correctness: " << (resid<SIM_THRES) << endl;
    return 0;
}

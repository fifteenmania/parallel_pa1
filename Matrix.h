#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <iomanip>
#include <pthread.h>
#include <cstring>
#include <math.h>
#include <limits>
#include <cstddef>
#define MAX_THREADS 512
#define SIM_THRES 0.0001
#define BLK_SIZE 32

class Matrix;

class Vector{
private:
    int n;
    double *data;
public:
    Vector();
    Vector(double *, int);
    Vector(const Vector &V);
    ~Vector();
    Vector &operator = (const Vector &other);
    // initialization
    int init_zeros(int);
    int init_ones(int);
    int init_rand(int, int);
    inline double *read_data();
    bool is_equal(const Vector &V);
    double residual_diff(const Vector &V);
    bool is_similar(const Vector &V);
    void print();
    void swaprow(int, int);
    void multrow(int, double);
    void submultrow(int, int, double);
    // Vector*Matrix
    void set_vmult(Matrix *, Vector *);
    Vector vmult(Matrix *);
};

class Matrix{
private:
    int n;
    double *data;
public:
    // static variables for multithread multiplication
    static double *Ad;
    static double *Bd;
    static double *Cd;
    static int size;
    static int num_threads;
    
    Matrix();
    Matrix(double *data, int n);
    Matrix(const Matrix &M);
    ~Matrix();
    Matrix &operator = (const Matrix &other);
    // initialization
    int init_zeros(int);
    int init_ones(int);
    int init_rand(int, int);
    inline double *read_data();
    void print();
    bool is_equal(const Matrix &B);
    bool is_similar(const Matrix &B);
    // matrix multiplication (with various algorithms)
    void set_mmult_naive(Matrix &A, Matrix &B);
    static void *_set_mmult_par_help(void *);
    void set_mmult_par(Matrix &A, Matrix &B, int p);
    void set_mmult_ordered(Matrix &A, Matrix &B);
    void set_mmult_blk(Matrix &A, Matrix &B);
    Matrix mmult(Matrix *A, Matrix *B);
    // Gaussian elimination
    void swaprow(int, int);
    void multrow(int, double);
    void submultrow(int, int, double);
    int pivotrow(int);
    void set_gauss_elim(Vector *);
    void set_backsub(Vector *);
};


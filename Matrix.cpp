#include "Matrix.h"

using namespace std;


// ********************************************
// Vector class definition
// ***********************************************


Vector::Vector()
{
    this->n = 0;
    this->data = NULL;
}

Vector::Vector(double *data, int n)
{
    this->n = n;
    this->data = data;
}
Vector::Vector(const Vector &V)
{
    int n = V.n;
    double *Vd = V.data;
    double *data = (double *)malloc(n*sizeof(double));
    memcpy(data, Vd, n*sizeof(double));
    this->n = n;
    this->data = data;
}
/*
Vector::Vector(const Vector &V)
{
    int n = V.n;
    double *Vd = V.data;
    double *data = (double *)malloc(n*sizeof(double));
    for (int i=0; i<n; i++){
        data[i] = Vd[i];
    }
    this->n = n;
    this->data = data;
}
*/
Vector::~Vector()
{
    if (this->data){
        free(this->data);
    }
    this->data = NULL;
    this->n = 0;
}

Vector& Vector::operator = (const Vector& other)
{
    if (this != &other){
        int n = other.n;
        double *Vd = other.data;
        double *data = (double *)malloc(n*sizeof(double));
        memcpy(data, Vd, n*sizeof(double));
        this->n = n;
        this->data = data;
        return *this;
    }
    return *this;
}

int Vector::init_zeros(int n)
{
    this->n = n;
    if (this->data){
        free(this->data);
    }
    double *A = (double*) malloc(n*sizeof(double));
    if (!A){
        return -1;
    }
    for (int i=0; i<n; i++){
        A[i] = 0;
    }
    this->data = A;
    return 0;
}

int Vector::init_ones(int n)
{
    this->n = n;
    if (this->data){
        free(this->data);
    }
    double *A = (double*) malloc(n*sizeof(double));
    if (!A){
        return -1;
    }
    for (int i=0; i<n; i++){
        A[i] = (double)1;
    }
    this->data = A;
    return 0;
}


int Vector::init_rand(int n, int seed)
{
    this->n = n;
    if (this->data){
        free(this->data);
    }
    double *A = (double*) malloc(n*sizeof(double));
    if (!A){
        return -1;
    }
    srand48(seed);
    for (int i=0; i<n; i++){
        A[i] = drand48();
    }
    this->data = A;
    return 0;
}

double *Vector::read_data()
{
    return this->data;
}

double Vector::residual_diff(const Vector &V)
{
    int n = V.n;
    if (this->n != n){
        return std::numeric_limits<double>::max();
    }
    double residual;
    double residual_sum = 0;
    for (int i=0; i<n; i++){
        residual = V.data[i] - this->data[i];
        residual_sum += residual*residual;
    }
    return sqrt(residual_sum/n);
}

bool Vector::is_similar(const Vector &V)
{
    return residual_diff(V) < SIM_THRES;
}

void Vector::print()
{
    int n = this->n;
    cout << setprecision(3) << fixed;
    for (int i=0; i<n; i++){
        cout << data[i] << endl;
    }
    return;
}

void Vector::swaprow(int row1, int row2)
{
    double *V = this->data;
    double tmp;
    tmp = V[row1];
    V[row1] = V[row2];
    V[row2] = tmp;
    return;
}

void Vector::submultrow(int row1, int row2, double a)
{
    double *V = this->data;
    V[row2] = V[row2] - a*V[row1];
    return;
}

void Vector::set_vmult(Matrix *A, Vector *X)
{
    int n = this->n;
    double *Ad = A->read_data();
    double *Xd = X->data;
    double *Zd = this->data;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            Zd[i] += Ad[i*n+j]*Xd[j];
        }
    }
    return;
}

Vector Vector::vmult(Matrix *A)
{
    Vector result = Vector();
    result.init_zeros(this->n);
    result.set_vmult(A, this);
    return result;
}


// **********************************************
// Matrix class definition
// **************************************************


double *Matrix::Ad = NULL;
double *Matrix::Bd = NULL;
double *Matrix::Cd = NULL;
int Matrix::num_threads = 0;
int Matrix::size = 0;
int Matrix::ppvt = 0;
int Matrix::cur_col = 0;
pthread_barrier_t Matrix::barrier;

Matrix::Matrix() 
{
    this->n = 0;
    this->data = NULL;
}

Matrix::Matrix(double *data, int n)
{
    this->n = n;
    this->data = data;
}

Matrix::Matrix(const Matrix &M)
{
    int n = M.n;
    double *A;
    this->n = n;
    A = (double*) malloc(n*n*sizeof(double));
    memcpy(A, M.data, n*n*sizeof(double));
    this->data = A;
}

Matrix::~Matrix() 
{
    if (this->data){
        free(this->data);
    }
    this->data = NULL;
    this->n = 0;
}

Matrix &Matrix::operator = (const Matrix &other)
{
    if (this != &other){
        int n = other.n;
        double *A;
        this->n = n;
        A = (double*) malloc(n*n*sizeof(double));
        memcpy(A, other.data, n*n*sizeof(double));
        this->data = A;
        return *this;
    }
    return *this;
}

int Matrix::init_zeros(int n)
{
    this->n = n;
    if (this->data){
        free(this->data);
    }
    double *A = (double*) malloc(n*n*sizeof(double));
    if (!A){
        return -1;
    }
    for (int i=0; i<n*n; i++){
        A[i] = 0;
    }
    this->data = A;
    return 0;
}

int Matrix::init_ones(int n)
{
    this->n = n;
    if (this->data){
        free(this->data);
    }
    double *A = (double*) malloc(n*n*sizeof(double));
    if (!A){
        return -1;
    }
    for (int i=0; i<n*n; i++){
        A[i] = (double)1;
    }
    this->data = A;
    return 0;
}


int Matrix::init_rand(int n, int seed)
{
    this->n = n;
    if (this->data){
        free(this->data);
    }
    double *A = (double*) malloc(n*n*sizeof(double));
    if (!A){
        return -1;
    }
    srand48(seed);
    for (int i=0; i<n*n; i++){
        A[i] = drand48();
    }
    this->data = A;
    return 0;
}

double *Matrix::read_data()
{
    return this->data;
}

void Matrix::print()
{
    int n = this->n;
    double *A = this->data;
    cout << setprecision(10)<< fixed;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            cout << A[n*i+j] << " ";
        }
        cout << endl;
    }
    return;
}

bool Matrix::is_similar(const Matrix &B)
{
    int n = this->n;
    double *lhsd = this->data;
    double *rhsd = B.data;
    if (n != B.n){
        return false;
    }
    double residual_sum = 0, residual;
    for (int i=0; i<n*n; i++){
        residual = rhsd[i]-lhsd[i];
        residual_sum += residual*residual;
    }
    return (residual_sum/n/n) < SIM_THRES; 
}


bool Matrix::is_equal(const Matrix &B)
{
    int n = this->n;
    double *lhsd = this->data;
    double *rhsd = B.data;
    if (n != B.n){
        return false;
    }
    for (int i=0; i<n*n; i++){
        if (rhsd[i] != lhsd[i]){
            return false;
        }
    }
    return true;
}

void Matrix::set_mmult_naive(Matrix &A, Matrix &B)
{
    int n = this->n;
    double *Cd = this->data;
    double *Ad = A.data;
    double *Bd = B.data;
    if ((n!=A.n) || (n!=B.n)){
        return;
    }
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                Cd[i*n+j] += Ad[i*n +k]*Bd[n*k+j];
            }
        }
    }
    return;
}

void *Matrix::_set_mmult_par_help(void *tnum_p)
{
    double *A = Matrix::Ad;
    double *B = Matrix::Bd;
    double *C = Matrix::Cd;
    int nt = Matrix::num_threads;
    int n = Matrix::size;
    int tnum = *((int *)tnum_p);
    int st = (int)((double)(tnum*n)/(double)nt);
    int ed = (int)((double)((tnum+1)*n)/(double)nt);
    for (int i=st; i<ed; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                C[i*n+k] += A[i*n+j]*B[n*j+k];
            }
        }
    }
    return NULL;
}

/*void *Matrix::_set_mmult_naive_par_help(void *context_p)
{
    mmult_binder *context = (mmult_binder *)context_p;
    return (context->self)->_set_mmult_naive_par(context->A, context->B, context->tid);
}

void *Matrix::_set_mmult_naive_par(double *A, double *B, int tid)
{
    double *C = this->data;
    int n = this->n;
    int st = (int)((double)(tid*n)/(double)NUM_THREADS);
    int ed = (int)((double)((tid+1)*n)/(double)NUM_THREADS);
    for (int i=st; i<ed; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                C[i*n+j] += A[i*n+k]*B[n*k+j];
            }
        }
    }
    return NULL;
}*/

void Matrix::set_mmult_par(Matrix &A, Matrix &B, int p)
{
    pthread_t p_threads[MAX_THREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    int tnum[MAX_THREADS];
    if ((n!=A.n) || (n!=B.n)){
        return;
    }
    Matrix::Ad = A.data;
    Matrix::Bd = B.data;
    Matrix::Cd = this->data;
    Matrix::num_threads = p;
    Matrix::size = this->n;
    for (int i=0; i<p; i++){
        tnum[i] = i;
        pthread_create(&p_threads[i], &attr, &Matrix::_set_mmult_par_help, (void *)&tnum[i]);
    }
    for (int j=0; j<p; j++){
        pthread_join(p_threads[j], NULL);
    }
    Matrix::Ad = NULL;
    Matrix::Bd = NULL;
    Matrix::Cd = NULL;
    Matrix::num_threads = 0;
    Matrix::size = 0;
    return;
}

void Matrix::set_mmult_ordered(Matrix &A, Matrix &B)
{
    int n = A.n;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            for (int k=0; k<n; k++){
                this->data[i*n+k] += A.data[i*n+j]*B.data[n*j+k];
            }
        }
    }
}

void Matrix::set_mmult_blk(Matrix &A, Matrix &B)
{
    // only be used when n is multiplier of 32
    int n = A.n;
    for (int i=0; i<n; i+=BLK_SIZE){
        for (int j=0; j<n; j+=BLK_SIZE){
            for (int k=0; k<n; k+=BLK_SIZE){
                // block loop
                for (int ii=i; ii<(i+BLK_SIZE); ii++){
                    for (int jj=j; jj<(j+BLK_SIZE); jj++){
                        for (int kk=k; kk<(k+BLK_SIZE); kk++){
                            this->data[ii*n+kk] += A.data[ii*n+jj]*B.data[n*jj+kk];
                        }
                    }
                }
            }
        }
    }
}

// assistant functions for Gaussian elimination
void Matrix::swaprow(int row1, int row2)
{
    // swap row1 and row2
    int n = this->n;
    double *A = this->data;
    double tmp;
    int r1_base = row1*n;
    int r2_base = row2*n;
    for (int i=0; i<n; i++){
        tmp = A[r1_base+i];
        A[r1_base+i] = A[r2_base+i];
        A[r2_base+i] = tmp;
    }
    return;
}

void Matrix::multrow(int row, double a)
{
    // multiply row1 by a
    int n = this->n;
    double *A = this->data;
    int r_base = row*n;
    for (int i=0; i<n; i++){
        A[r_base+i] = A[r_base+i]*a;
    }
    return;
}


void Matrix::submultrow(int row1, int row2, double a)
{
    // subtract row2 by row1*a
    int n = this->n;
    double *A = this->data;
    int r1_base = row1*n;
    int r2_base = row2*n;
    for (int i=row1; i<n; i++){
        A[r2_base+i] = A[r2_base+i] - a*A[r1_base+i];
    }
    return;
}

int Matrix::pivotrow(int col)
{
    // find maximal row with given column(partial pivot)
    // return -1 if no nonzero entry.
    int n = this->n;
    double *A = this->data;
    int pivot_row = -1;
    double pivot_max = 0;
    for (int i=col; i<n; i++){
        if (fabs(A[i*n+col])>pivot_max){
            pivot_max = abs(A[i*n+col]);
            pivot_row = i;
        }
    }
    return pivot_row;
}

void *Matrix::_set_gauss_elim_par_help(void *tnum_p)
{
    double a;
    int n = size;
    int tnum = *((int*)tnum_p);
    for (int i=cur_col+tnum+1; i<n; i+=num_threads){
        a = Ad[i*n+cur_col]/Ad[cur_col*n+cur_col];
        for (int j=i; j<n; j++){
            Ad[i*n+j] -= a*Ad[cur_col*n+j];
        }
        Bd[i] -= a*Bd[cur_col];
    }
    return NULL;
}

void Matrix::set_gauss_elim_par(Vector *V, int p)
{
    // set this, V as echelon form
    // Modify V and this
    pthread_t p_threads[MAX_THREADS];
    int pivot_row;
    int n = this->n;
    int tnum[MAX_THREADS];
    for (int i=0; i<p; i++){
        tnum[i] = i;
    }
    Ad = this->data;
    Bd = V->data;
    size = n;
    num_threads = p;
    for (int i=0; i<n-1; i++){
        pivot_row = this->pivotrow(i);
        if (pivot_row != -1){
            ppvt = pivot_row;
            cur_col = i;
            this->swaprow(i, pivot_row);
            V->swaprow(i, pivot_row);
            for (int j=0; j<p; j++){
                pthread_create(&p_threads[j], NULL, Matrix::_set_gauss_elim_par_help, (void *)&(tnum[j]));
            }
            for (int j=0; j<p; j++){
                pthread_join(p_threads[j], NULL);
            }
            /*this->swaprow(i, pivot_row);
            V->swaprow(i, pivot_row);
            for (int j=i+1; j<n; j++){
                a = A[j*n+i]/A[i*n+i];
                this->submultrow(i,j,a);
                V->submultrow(i,j,a);
            }*/
        }
    }
    return;
}

void Matrix::set_gauss_elim(Vector *V)
{
    // set this, V as echelon form
    // Modify V and this
    int pivot_row;
    int n = this->n;
    double *A = this->data;
    double a;
    for (int i=0; i<n-1; i++){
        pivot_row = this->pivotrow(i);
        if (pivot_row != -1){
            this->swaprow(i, pivot_row);
            V->swaprow(i, pivot_row);
            for (int j=i+1; j<n; j++){
                a = A[j*n+i]/A[i*n+i];
                this->submultrow(i,j,a);
                V->submultrow(i,j,a);
            }
        }
    }
    return;
}

void *Matrix::_set_backsub_par_help(void *tnum_p)
{
    int n = Matrix::size;
    double *A = Matrix::Ad;
    double *X = Matrix::Bd;
    int p = Matrix::num_threads;
    int tnum = *(int*)tnum_p;
    double a;
    int st_idx, ed_idx;
    for (int i=n-1; i>=LINE_SIZE; i--){
        if (A[i*n+i] != 0){
            a = X[i]/A[i*n+i];
            X[i] = a;
            st_idx = (tnum*i)/p;
            ed_idx = ((tnum+1)*i)/p;
            for (int j=st_idx; j<ed_idx; j++){
                X[j] -= A[j*n+i]*a;
            }
        }
        pthread_barrier_wait(&(Matrix::barrier));
    } 
    return NULL;
}

void Matrix::set_backsub_par(Vector *V, int p)
{
    Matrix::Ad = this->data;
    Matrix::Bd = V->data;
    Matrix::size = this->n;
    Matrix::num_threads = p;

    int tnum[MAX_THREADS];
    for (int i=0; i<p; i++){
        tnum[i] = i;
    }
    
    pthread_t p_threads[MAX_THREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_barrier_init(&Matrix::barrier, NULL, p);
    for (int k=0; k<p; k++){
        pthread_create(&p_threads[k], &attr, &Matrix::_set_backsub_par_help, (void *)&tnum[k]);
    }
    for (int l=0; l<p; l++){
        pthread_join(p_threads[l], NULL);
    }
    double a;
    // do remaining small tasks in single thread mode
    for (int m=min(LINE_SIZE, n-1); m>=0; m--){
        if (Ad[m*n+m] != 0){
            a = Bd[m]/Ad[m*n+m];
            Bd[m] = a;
            for (int q=m-1; q>=0; q--){
                Bd[q] -= Ad[q*n+m]*a;
            }
        }
    }

    Matrix::Ad = NULL;
    Matrix::Bd = NULL;
    Matrix::size = 0;
    Matrix::num_threads = 0;
    return;
}

void Matrix::set_backsub(Vector *V)
{
    // set V as X vector
    // modify V
    int n = this->n;
    double *A = this->data;
    double *Xd = V->read_data();
    double a;
    for (int i=n-1; i>=0; i--){
        if (A[i*n+i] != 0){
            a = Xd[i]/A[i*n+i];
            Xd[i] = a;
            for (int j=i-1; j>=0; j--){
                Xd[j] -= A[j*n+i]*a;
            }
        }
    }       
    return;
}



#ifndef matrix_mult_mt_h
#define matrix_mult_mt_h

#include <TMatrixD.h>

// data structure to pass to threads
struct input_data {

    const double* a_array_;
    const double* b_array_;
    double* c_array_;
    int a_nrow_;
    int a_ncol_;
    int b_nrow_;
    int b_ncol_;
    int start_;
    int end_;
        
};


// function that's passed to each thread
void* multiply_on_thread(void* data);


// matrix multiplication where matrix data are stored in arrays
void multiply_mt_impl(const double* a_array, const double* b_array, double* c_array,
                      const int a_nrow, const int a_ncol, const int b_nrow, const int b_ncol);

// matrix multiplication using root TMatrixD objects
TMatrixD multiply_mt(const TMatrixD& a, const TMatrixD& b);

// Perform matrix multiplication: C = A B A^T using multithread
void ABAT_MT(const TMatrixD& a, const TVectorD& b, TMatrixD& c);
  
#endif

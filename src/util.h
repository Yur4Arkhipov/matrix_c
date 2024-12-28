#ifndef S21_UTIL
#define S21_UTIL

#include "s21_matrix.h"

int s21_is_matrix_valid(matrix_t *A);
void shr_matrix(matrix_t *A, matrix_t *minor, int row_to_skip, int col_to_skip);
void get_minor_matrix(matrix_t *A, matrix_t *minor, matrix_t *result, int size);
void calc_determinant(matrix_t *A, double *determinant_sum, int col,
                      double minor_det);
#endif

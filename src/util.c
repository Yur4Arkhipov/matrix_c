#include "util.h"

#include "s21_matrix.h"

int s21_is_matrix_valid(matrix_t *A) {
  int return_value = 1;
  if (!A || A->rows <= 0 || A->columns <= 0 || !A->matrix) {
    return_value = 0;
  }

  if (return_value) {
    for (int i = 0; i < A->rows; i++) {
      if (!A->matrix[i]) {
        return_value = 0;
      }
    }
  }

  return return_value;
}

void shr_matrix(matrix_t *A, matrix_t *minor, int row_to_skip,
                int col_to_skip) {
  int minor_row = 0;
  for (int row = 0; row < A->rows; row++) {
    if (row == row_to_skip) continue;

    int minor_col = 0;
    for (int col = 0; col < A->columns; col++) {
      if (col == col_to_skip) continue;

      minor->matrix[minor_row][minor_col] = A->matrix[row][col];
      minor_col++;
    }
    minor_row++;
  }
}

void get_minor_matrix(matrix_t *A, matrix_t *minor, matrix_t *result,
                      int size) {
  int status;

  for (int row = 0; row < size; row++) {
    for (int col = 0; col < size; col++) {
      shr_matrix(A, minor, row, col);

      double minor_det = 0;
      status = s21_determinant(minor, &minor_det);
      if (status == OK) {
        if ((row + col) % 2 == 1 && minor_det != 0) {
          minor_det *= -1;
        }
        result->matrix[row][col] = minor_det;
      }
    }
  }
}

void calc_determinant(matrix_t *A, double *determinant_sum, int col,
                      double minor_det) {
  if (col % 2 == 0) {
    *determinant_sum += A->matrix[0][col] * minor_det;
  } else {
    *determinant_sum -= A->matrix[0][col] * minor_det;
  }
}
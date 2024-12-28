#include "s21_matrix.h"

#include "util.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = OK;

  if (!result || rows < 1 || columns < 1) {
    return INVALID_MATRIX;
  }

  result->rows = rows;
  result->columns = columns;

  double **matrix = malloc(rows * sizeof(double *));
  if (matrix) {
    int allocation_failed = 0;
    for (int i = 0; i < rows && !allocation_failed; i++) {
      matrix[i] = malloc(columns * sizeof(double));
      if (matrix[i] == NULL) {
        allocation_failed = 1;
        status = INVALID_MATRIX;

        for (int j = 0; j < i; j++) {
          free(matrix[j]);
          matrix[j] = NULL;
        }
        free(matrix);
        matrix = NULL;
        result->matrix = NULL;
        result->rows = 0;
        result->columns = 0;
      }
    }

    if (!allocation_failed) {
      result->matrix = matrix;
    }
  } else {
    status = INVALID_MATRIX;
  }

  return status;
}

void s21_remove_matrix(matrix_t *A) {
  if (!A || !A->matrix) {
    return;
  }

  for (int i = 0; i < A->rows; ++i) {
    if (A->matrix[i]) {
      free(A->matrix[i]);
    }
  }
  free(A->matrix);
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (!s21_is_matrix_valid(A) || !s21_is_matrix_valid(B)) {
    return FAILURE;
  }

  int status = SUCCESS;

  if (A->rows != B->rows || A->columns != B->columns) {
    status = FAILURE;
  } else {
    for (int i = 0; i < A->rows && status; i++) {
      for (int j = 0; j < A->columns && status; j++) {
        status = fabs(A->matrix[i][j] - B->matrix[i][j]) < 1e-7;
      }
    }
  }

  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;

  if (!result || !s21_is_matrix_valid(A) || !s21_is_matrix_valid(B)) {
    return INVALID_MATRIX;
  }

  if (A->rows != B->rows || A->columns != B->columns) {
    status = CALCULATION_ERROR;
  } else {
    status = s21_create_matrix(A->rows, A->columns, result);
  }

  if (status == OK) {
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
      }
    }
  }

  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (!result || !s21_is_matrix_valid(A) || !s21_is_matrix_valid(B)) {
    return INVALID_MATRIX;
  }

  int status = OK;

  if (A->rows != B->rows || A->columns != B->columns) {
    status = CALCULATION_ERROR;
  } else {
    status = s21_create_matrix(A->rows, A->columns, result);
  }

  if (status == OK) {
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
      }
    }
  }

  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (!result || !s21_is_matrix_valid(A)) {
    return INVALID_MATRIX;
  }

  int status = s21_create_matrix(A->rows, A->columns, result);

  if (status == OK) {
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[i][j] = number * A->matrix[i][j];
      }
    }
  }

  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;

  if (!result || !s21_is_matrix_valid(A) || !s21_is_matrix_valid(B)) {
    return INVALID_MATRIX;
  }

  if (A->columns != B->rows) {
    status = CALCULATION_ERROR;
  } else {
    status = s21_create_matrix(A->rows, B->columns, result);
  }
  if (status == OK) {
    for (int i = 0; i < result->rows; i++) {
      for (int j = 0; j < result->columns; j++) {
        result->matrix[i][j] = 0.0;

        for (int k = 0; k < A->columns; k++) {
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
        }
      }
    }
  }

  return status;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (!result || !s21_is_matrix_valid(A)) {
    return INVALID_MATRIX;
  }

  int status = s21_create_matrix(A->columns, A->rows, result);

  if (status == OK) {
    for (int i = 0; i < A->rows; ++i) {
      for (int j = 0; j < A->columns; ++j) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }

  return status;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (!result || !s21_is_matrix_valid(A)) {
    return INVALID_MATRIX;
  }
  int status = OK;
  int size = A->rows;

  if (A->rows != A->columns) {
    status = CALCULATION_ERROR;
    result->matrix = NULL;
  } else {
    status = s21_create_matrix(size, size, result);
  }

  if (status == OK) {
    if (size == 1) {
      result->matrix[0][0] = 1;
    } else {
      matrix_t minor;
      status = s21_create_matrix(size - 1, size - 1, &minor);
      if (status == OK) {
        get_minor_matrix(A, &minor, result, size);
        s21_remove_matrix(&minor);
      } else {
        s21_remove_matrix(result);
      }
    }
  }

  return status;
}

int s21_determinant(matrix_t *A, double *result) {
  if (!result || !s21_is_matrix_valid(A)) {
    return INVALID_MATRIX;
  }

  int status = OK;
  int size = A->rows;

  if (A->rows != A->columns) {
    status = CALCULATION_ERROR;
  } else if (size == 1) {
    *result = A->matrix[0][0];
  } else if (size == 2) {
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[1][0] * A->matrix[0][1];
  } else {
    matrix_t minor;
    double determinant_sum = 0;

    status = s21_create_matrix(size - 1, size - 1, &minor);
    if (status == OK) {
      for (int col = 0; col < size && status == OK; col++) {
        shr_matrix(A, &minor, 0, col);

        double minor_det = 0;
        status = s21_determinant(&minor, &minor_det);
        if (status == OK) {
          calc_determinant(A, &determinant_sum, col, minor_det);
        }
      }

      *result = determinant_sum;
    }

    s21_remove_matrix(&minor);
  }

  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  if (!result || !s21_is_matrix_valid(A)) {
    return INVALID_MATRIX;
  }

  int status = OK;
  int size = A->rows;

  if (A->rows != A->columns) {
    return CALCULATION_ERROR;
  }

  matrix_t complements;
  matrix_t transposed;
  complements.matrix = NULL;
  transposed.matrix = NULL;
  double det = 0;

  status = s21_determinant(A, &det);
  if (status == OK && fabs(det) < 1e-7) {
    return CALCULATION_ERROR;
  }

  if (status == OK) {
    status = s21_calc_complements(A, &complements);
  }

  if (status == OK) {
    status = s21_transpose(&complements, &transposed);
  }

  if (status == OK) {
    status = s21_create_matrix(size, size, result);
  }

  if (status == OK) {
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        result->matrix[i][j] = transposed.matrix[i][j] / det;
      }
    }
  }

  if (complements.matrix) {
    s21_remove_matrix(&complements);
  }

  if (transposed.matrix) {
    s21_remove_matrix(&transposed);
  }

  return status;
}

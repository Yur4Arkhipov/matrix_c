#ifndef TESTS_H
#define TESTS_H

#include "../s21_matrix.h"

#define TEST_SIZE 625

typedef struct {
  int rows;
  int cols;
  int start;
} test_matrix_t;

typedef struct {
  test_matrix_t matrix;
  int code;
} test_matrix_code_t;

typedef struct {
  test_matrix_t a;
  test_matrix_t b;

  int eq_result;

  test_matrix_code_t sum_result;
  test_matrix_code_t sub_result;

  double mult_number;
  test_matrix_t mult_number_result;
  test_matrix_code_t mult_matrix_result;

  test_matrix_t transpose_result;

  test_matrix_code_t calc_complements_result;
  int determinant_code;
  double determinant_result;
  test_matrix_code_t inverse_result;
} dataset_t;

extern const double dataset_values[];
extern const dataset_t dataset_tests[];

#endif

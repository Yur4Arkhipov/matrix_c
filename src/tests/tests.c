#include "tests.h"

#include <check.h>
#include <stdio.h>

matrix_t get_matrix(test_matrix_t A) {
  matrix_t a;
  s21_create_matrix(A.rows, A.cols, &a);
  int idx = A.start;
  for (int i = 0; i < A.rows; i++) {
    for (int j = 0; j < A.cols; j++) {
      a.matrix[i][j] = dataset_values[idx];
      idx++;
    }
  }
  return a;
}

int get_matrix_code(test_matrix_code_t A, matrix_t* result) {
  int code = A.code;
  if (!code) {
    *result = get_matrix(A.matrix);
  } else {
    result->matrix = NULL;
  }
  return code;
}

void run_test_eq(matrix_t* A, matrix_t* B, const dataset_t* test) {
  int result = s21_eq_matrix(A, B);
  ck_assert_int_eq(result, test->eq_result);
}

void run_test_sum(matrix_t* A, matrix_t* B, const dataset_t* test) {
  matrix_t result_test;
  int code_test = get_matrix_code(test->sum_result, &result_test);
  matrix_t result;
  int code = s21_sum_matrix(A, B, &result);
  ck_assert_int_eq(code_test, code);
  if (!code) {
    ck_assert(s21_eq_matrix(&result, &result_test));
    s21_remove_matrix(&result);
  }
  s21_remove_matrix(&result_test);
}

void run_test_sub(matrix_t* A, matrix_t* B, const dataset_t* test) {
  matrix_t result_test;
  int code_test = get_matrix_code(test->sub_result, &result_test);
  matrix_t result;
  int code = s21_sub_matrix(A, B, &result);
  ck_assert_int_eq(code_test, code);
  if (!code) {
    ck_assert(s21_eq_matrix(&result, &result_test));
    s21_remove_matrix(&result);
  }
  s21_remove_matrix(&result_test);
}

void run_test_mult_number(matrix_t* A, const dataset_t* test) {
  matrix_t result_test = get_matrix(test->mult_number_result);
  matrix_t result;
  int code = s21_mult_number(A, test->mult_number, &result);
  ck_assert_int_eq(code, OK);
  ck_assert(s21_eq_matrix(&result, &result_test));
  s21_remove_matrix(&result_test);
  s21_remove_matrix(&result);
}

void run_test_mult_matrix(matrix_t* A, matrix_t* B, const dataset_t* test) {
  matrix_t result_test;
  int code_test = get_matrix_code(test->mult_matrix_result, &result_test);
  matrix_t result;
  int code = s21_mult_matrix(A, B, &result);
  ck_assert_int_eq(code_test, code);
  if (!code) {
    ck_assert(s21_eq_matrix(&result, &result_test));
    s21_remove_matrix(&result);
  }
  s21_remove_matrix(&result_test);
}

void run_test_transpose(matrix_t* A, const dataset_t* test) {
  matrix_t result_test = get_matrix(test->transpose_result);
  matrix_t result;
  int code = s21_transpose(A, &result);
  ck_assert_int_eq(code, OK);
  ck_assert(s21_eq_matrix(&result, &result_test));
  s21_remove_matrix(&result_test);
  s21_remove_matrix(&result);
}

void run_test_calc_complements(matrix_t* A, const dataset_t* test) {
  matrix_t result_test;
  int code_test = get_matrix_code(test->calc_complements_result, &result_test);
  matrix_t result;
  int code = s21_calc_complements(A, &result);
  ck_assert_int_eq(code_test, code);
  if (!code) {
    ck_assert(s21_eq_matrix(&result, &result_test));
    s21_remove_matrix(&result);
  }
  s21_remove_matrix(&result_test);
}

void run_test_determinant(matrix_t* A, const dataset_t* test) {
  double result_test = test->determinant_result;
  int code_test = test->determinant_code;
  double result;
  int code = s21_determinant(A, &result);
  ck_assert_int_eq(code_test, code);
  if (!code) {
    ck_assert_double_eq_tol(result, result_test, 1e-7);
  }
}

void run_test_inverse(matrix_t* A, const dataset_t* test) {
  matrix_t result_test;
  int code_test = get_matrix_code(test->inverse_result, &result_test);
  matrix_t result;
  int code = s21_inverse_matrix(A, &result);
  ck_assert_int_eq(code_test, code);
  if (!code) {
    ck_assert(s21_eq_matrix(&result, &result_test));
    s21_remove_matrix(&result);
  }
  s21_remove_matrix(&result_test);
}

void run_test(int idx) {
  const dataset_t* test = &dataset_tests[idx];

  matrix_t A = get_matrix(test->a);
  matrix_t B = get_matrix(test->b);

  run_test_eq(&A, &B, test);
  run_test_sum(&A, &B, test);
  run_test_sub(&A, &B, test);
  run_test_mult_number(&A, test);
  run_test_mult_matrix(&A, &B, test);
  run_test_transpose(&A, test);
  run_test_calc_complements(&A, test);
  run_test_determinant(&A, test);
  run_test_inverse(&A, test);

  s21_remove_matrix(&A);
  s21_remove_matrix(&B);
}

START_TEST(test_all_valid) { run_test(_i); }
END_TEST

START_TEST(test_all_invalid) {
  ck_assert_int_eq(s21_eq_matrix(NULL, NULL), 0);
  ck_assert_int_eq(s21_sum_matrix(NULL, NULL, NULL), INVALID_MATRIX);
  ck_assert_int_eq(s21_sub_matrix(NULL, NULL, NULL), INVALID_MATRIX);
  ck_assert_int_eq(s21_mult_number(NULL, 1.0, NULL), INVALID_MATRIX);
  ck_assert_int_eq(s21_mult_matrix(NULL, NULL, NULL), INVALID_MATRIX);
  ck_assert_int_eq(s21_transpose(NULL, NULL), INVALID_MATRIX);
  ck_assert_int_eq(s21_calc_complements(NULL, NULL), INVALID_MATRIX);
  ck_assert_int_eq(s21_determinant(NULL, NULL), INVALID_MATRIX);
  ck_assert_int_eq(s21_inverse_matrix(NULL, NULL), INVALID_MATRIX);

  matrix_t m;
  s21_create_matrix(3, 3, &m);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      m.matrix[i][j] = 0.0;
    }
  }
  matrix_t result;
  ck_assert_int_eq(s21_inverse_matrix(&m, &result), CALCULATION_ERROR);
  s21_remove_matrix(&m);
}
END_TEST

int main() {
  Suite* s = suite_create("test_all");
  TCase* tc = tcase_create("test_all");
  tcase_add_loop_test(tc, test_all_valid, 0, TEST_SIZE);
  tcase_add_test(tc, test_all_invalid);

  suite_add_tcase(s, tc);
  SRunner* sr = srunner_create(s);
  srunner_set_fork_status(sr, CK_NOFORK);
  srunner_run_all(sr, CK_NORMAL);
  srunner_free(sr);
}
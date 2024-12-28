CC=gcc
CC_FLAGS=-c -std=c11 -Wall -Werror -Wextra
GCOV=-fprofile-arcs -ftest-coverage
TEST_FLAGS=-lcheck -lsubunit -lrt -lm -lpthread

# clang-format
CLANG_FORMAT = clang-format -i tests/*.c tests/*.h *.c *.h
CLANG_FORMAT_CHECK = clang-format -n tests/*.c tests/*.h *.c *.h

all: build test gcov_report

build: s21_matrix.a

test: s21_matrix.c util.c tests/generated.c tests/tests.c
	${CC} s21_matrix.c util.c tests/generated.c tests/tests.c ${TEST_FLAGS} -o tests/test
	./tests/test

rebuild: clean build

s21_matrix.a: s21_matrix.o util.o
	ar rcs s21_matrix.a s21_matrix.o util.o
	ranlib s21_matrix.a

s21_matrix.o: s21_matrix.c
	${CC} ${CC_FLAGS} s21_matrix.c

s21_util.o: util.c
	${CC} ${CC_FLAGS} util.c

gcov_report: s21_matrix.c util.c tests/generated.c tests/tests.c
	${CC} --coverage s21_matrix.c util.c tests/generated.c tests/tests.c ${TEST_FLAGS} -o tests/test_report
	./tests/test_report
	lcov --directory . --capture -o coverage.info
	genhtml --output-directory report --legend coverage.info
	rm ./tests/*.gcda
	rm ./tests/*.gcno
	rm coverage.info

clang_format: 
	${CLANG_FORMAT}
	${CLANG_FORMAT_CHECK}

clean:
	-rm -rf ./*.o
	-rm ./*.a
	-rm ./tests/test
	-rm ./tests/test_report
	-rm -rf report
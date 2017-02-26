CC=g++


serial:
	$(CC) -o main_serial main_serial.cpp
	./main_serial

utest:
	$(CC) -o main_unit_test main_unit_test.cpp
	./main_unit_test

vtest:
	$(CC) -o main_verification_test main_verification_test.cpp -fopenmp
	./main_verification_test

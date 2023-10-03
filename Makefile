CC=icpc

WARNINGS = -Wall -Wextra -Wunused -Woverloaded-virtual -Wshadow \
           -Wpointer-arith -Wcast-qual -Wwrite-strings -Wconversion \
           -wd111,304,383,967,981,1418,1419,1572,1682,2259 \
           -diag-disable ipo -diag-disable=10441

CFLAGS= $(WARNINGS) -O3 -xHost -fp-model precise -fp-model except -fp-model source -std=c++11 -qmkl=sequential

run: TestLinearSolvers
		echo "Running with tolerance 1E-15"
		./TestLinearSolvers

lower_tolerance:
		sed -i "s/const double tolerance = 1E-15;/const double tolerance = 1E-13;/g" ./TestLinearSolvers.cpp
reset_tolerance:
		sed -i "s/const double tolerance = 1E-13;/const double tolerance = 1E-15;/g" ./TestLinearSolvers.cpp 

run_lower_tol: lower_tolerance all reset_tolerance
		echo "Running with tolerance 1E-13"
		./TestLinearSolvers

all: TestLinearSolvers

OBJECTS = TestLinearSolvers.o LinearSolvers.o

TestLinearSolvers: $(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) -o TestLinearSolvers

%.o: %.cpp
		$(CC) $(CFLAGS) -c $<

clean:
		rm -rf *.o
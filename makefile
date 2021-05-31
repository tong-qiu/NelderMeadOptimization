output: NumericalOptimization.o main.o
	g++ -std=c++11 -g NumericalOptimization.o main.o -o output
NumericalOptimization.o: NumericalOptimization.cpp NumericalOptimization.h
	g++ -std=c++11 -g -c NumericalOptimization.cpp
main.o: main.cpp     
	g++ -std=c++11 -g -c main.cpp
clean:
	rm NumericalOptimization.o main.o

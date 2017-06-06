CFLAGS=`pkg-config --cflags gsl` -g -Wall -O0
LOADLIBES=`pkg-config --libs gsl` -lm

src/mixture_bernoulli: src/metropolis.o src/random.o  src/common.o
src/multinomial_example: src/metropolis.o src/random.o  src/gmm.o src/discrete.o src/common.o src/gsl.o 
src/normal_example: src/metropolis.o src/random.o  src/gmm.o src/common.o src/gsl.o
clean:
	rm src/*.o

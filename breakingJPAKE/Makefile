all:
	g++ \
		-g -O2 -std=c++11 -pthread -march=native \
		jpake.c \
		bench.c \
		-o main \
		-lntl -lgmp -lm 

clean:
	@rm -f main bench


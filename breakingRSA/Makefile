all:
	g++ \
		-g -O2 -std=c++11 -pthread -march=native \
		rsa.c \
		bench.c \
		-o main \
		-lntl -lgmp -lm 

clean:
	@rm -f main bench


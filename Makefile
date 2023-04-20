CC = g++
CFLAGS = -Wall -O3 -std=c++17 -fopenmp
OUT = ising
graph ?= ising
build:
	$(CC) $(CFLAGS) main.cpp ising.cpp -o $(OUT)

run:
	./$(OUT)

plot:
	python plot.py -f out/$(graph).csv -o graphs/$(graph).png
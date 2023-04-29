CC = g++
CFLAGS = -Wall -O3 -std=c++17 -fopenmp
OUT = ising
name ?= ising
build:
	$(CC) $(CFLAGS) main.cpp ising.cpp -o $(OUT)

run:
	./$(OUT)

plot:
	python plot.py -f out/$(name).csv -o graphs/$(name).png
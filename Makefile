all: build
	g++-16 ode/main.cpp -o build/main -std=c++23

build:
	mkdir build

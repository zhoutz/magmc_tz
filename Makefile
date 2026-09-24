all: build
	g++-16 src/bench/base.cpp -o build/base -std=c++23

build:
	mkdir build

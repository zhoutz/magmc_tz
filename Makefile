all: build
	g++-16 src/bench/base.cpp -o build/base -std=c++23
	g++-16 src/bench/ref.cpp -o build/ref -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib

build:
	mkdir build

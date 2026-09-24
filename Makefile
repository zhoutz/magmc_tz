# all: build
# 	g++-16 src/bench/base.cpp -o build/base -std=c++23
# 	g++-16 src/bench/ref.cpp -o build/ref -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib
	
base:
	g++-16 src/bench/base.cpp -o build/base -std=c++23 -O3
	build/base
	
ref:
	g++-16 src/bench/ref.cpp -o build/ref -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib -O3
	build/ref

build:
	mkdir build

# all: build
# 	g++-16 src/bench/base.cpp -o build/base -std=c++23
# 	g++-16 src/bench/ref.cpp -o build/ref -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib
	
base:
	g++-16 src/bench/base.cpp -o build/base -std=c++23 -O3
	build/base
	
ref:
	g++-16 src/bench/ref.cpp -o build/ref -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib -O3
	time build/ref
	
human:
	g++-16 src/bench/human.cpp -o build/human -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib -O3
	time build/human
	
main:
	g++-16 src/dev/main.cpp -o build/main -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib -O3
	time build/main
	
ft07:
	g++-16 src/dev/ft07.cpp -o build/ft07 -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib -O3
	time build/ft07 --photons 1000000

build:
	mkdir build

.PHONY: fd11 fd11-test
fd11:
	mkdir -p build output
	g++-16 src/dev/fd11.cpp -o build/fd11 -std=c++23 -O3 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib

fd11-test: fd11
	g++-16 src/dev/fd11_test.cpp -o build/fd11_test -std=c++23 -O3 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib
	build/fd11_test

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

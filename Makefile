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
	
	
test_qagp:
	g++-16 src/bench/test_qagp.cpp -o build/test_qagp -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib -O3
	time build/test_qagp

build:
	mkdir build

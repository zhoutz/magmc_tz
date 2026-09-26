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

.PHONY: global_sqrt velocity_panels transport_fast benchmark_optical_depth
global_sqrt velocity_panels transport_fast:
	mkdir -p build output
	g++-16 src/bench/$@.cpp -o build/$@ -std=c++23 -O3 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib
	build/$@

benchmark_optical_depth:
	python3 py/compare_all_optical_depth.py

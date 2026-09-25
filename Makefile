# all: build
# 	g++-16 src/bench/base.cpp -o build/base -std=c++23
# 	g++-16 src/bench/ref.cpp -o build/ref -std=c++23 -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib

CXXFLAGS = -std=c++23 -O3
GSL_FLAGS = -lgsl -I/opt/homebrew/include -L/opt/homebrew/lib
CXX = g++-16

base:
	$(CXX) src/bench/base.cpp -o build/base $(CXXFLAGS)
	build/base

ref:
	$(CXX) src/bench/ref.cpp -o build/ref $(CXXFLAGS) $(GSL_FLAGS)
	build/ref

step_limiter:
	$(CXX) src/bench/step_limiter.cpp -o build/step_limiter $(CXXFLAGS) $(GSL_FLAGS)
	build/step_limiter

gauss_legendre:
	$(CXX) src/bench/gauss_legendre.cpp -o build/gauss_legendre $(CXXFLAGS)
	build/gauss_legendre

gauss_legendre_10:
	$(CXX) src/bench/gauss_legendre_10.cpp -o build/gauss_legendre_10 $(CXXFLAGS)
	build/gauss_legendre_10

adaptive_gauss:
	$(CXX) src/bench/adaptive_gauss.cpp -o build/adaptive_gauss $(CXXFLAGS)
	build/adaptive_gauss

event_bracketed:
	$(CXX) src/bench/event_bracketed.cpp -o build/event_bracketed $(CXXFLAGS) $(GSL_FLAGS)
	build/event_bracketed

test_nonradial:
	$(CXX) src/bench/test_nonradial.cpp -o build/test_nonradial $(CXXFLAGS) $(GSL_FLAGS)
	build/test_nonradial

all_methods: step_limiter gauss_legendre gauss_legendre_10 adaptive_gauss

compare: all_methods
	@echo "Comparing all methods against benchmark..."
	@python3 py/cmp_od.py step_limiter
	@python3 py/cmp_od.py gauss_legendre
	@python3 py/cmp_od.py gauss_legendre_10
	@python3 py/cmp_od.py adaptive_gauss

build:
	mkdir build

clean:
	rm -rf build/*.o build/base build/ref build/step_limiter build/gauss_legendre* build/adaptive_gauss build/event_bracketed build/test_nonradial

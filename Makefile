ifeq ($(origin CXX),default)
CXX := g++-16
endif
CXXFLAGS ?= -O3 -std=c++23 -Wall -Wextra -Wpedantic
HEADERS := $(wildcard ode/*.hpp)
TESTS := test_dopr5 reference_probe test_null_geodesic test_transport_smoke test_fast_transport

.PHONY: all check benchmark benchmark-quick benchmark-fast benchmark-ensemble
all: build/main

build:
	mkdir -p build

build/main: ode/main.cpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) $< -o $@

build/benchmark_resonance: test/benchmark_resonance.cpp test/benchmark_cases.hpp test/resonance_reference.hpp $(HEADERS) | build
	$(CXX) $(CXXFLAGS) $< -o $@

build/%: test/%.cpp $(wildcard test/*.hpp) $(HEADERS) | build
	$(CXX) $(CXXFLAGS) -Iode $< -o $@

check: $(addprefix build/,$(TESTS))
	./build/test_dopr5
	./build/reference_probe
	./build/test_null_geodesic
	./build/test_transport_smoke
	./build/test_fast_transport

benchmark: build/benchmark_resonance
	./build/benchmark_resonance --repeats 5 --output benchmark/results.csv
	python3 benchmark/summarize.py benchmark/results.csv

benchmark-quick: build/benchmark_resonance
	./build/benchmark_resonance --quick --repeats 3 --output benchmark/quick.csv
	python3 benchmark/summarize.py benchmark/quick.csv

benchmark-fast: build/benchmark_resonance
	./build/benchmark_resonance --fast-comparison --repeats 7 --output benchmark/fast_results.csv
	python3 benchmark/summarize_fast.py

benchmark-ensemble: build/main
	python3 benchmark/benchmark_ensemble.py

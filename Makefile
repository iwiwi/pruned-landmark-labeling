CXX = g++
CXXFLAGS = -g -Wall -Wextra -O3


all: bin bin/construct_index bin/query_distance bin/benchmark

bin:
	mkdir -p bin

bin/construct_index: samples/construct_index_main.cc src/pruned_landmark_labeling.h
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin/query_distance: samples/query_distance_main.cc src/pruned_landmark_labeling.h
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin/benchmark: samples/benchmark_main.cc src/pruned_landmark_labeling.h
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin/pruned_landmark_labeling_test: src/pruned_landmark_labeling_test.cc src/pruned_landmark_labeling.h
	$(CXX) $(CXXFLAGS) -lgtest -lgtest_main -o $@ $^



.PHONY: test clean

test: bin/pruned_landmark_labeling_test
	bin/pruned_landmark_labeling_test

clean:
	rm -rf bin
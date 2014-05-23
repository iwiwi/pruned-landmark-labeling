#include <sys/time.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "pruned_landmark_labeling.h"

using namespace std;

const int kNumQueries = 1000000;

double GetCurrentTimeSec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: construct_index GRAPH" << endl;
  }

  PrunedLandmarkLabeling<0> pll;
  if (!pll.ConstructIndex(argv[1])) {
    cerr << "error: Load or construction failed" << endl;
    exit(EXIT_FAILURE);
  }
  pll.PrintStatistics();

  vector<int> vs(kNumQueries), ws(kNumQueries);
  for (int i = 0; i < kNumQueries; ++i) {
    vs[i] = rand() % pll.GetNumVertices();
    ws[i] = rand() % pll.GetNumVertices();
  }

  double time_start = GetCurrentTimeSec();
  for (int i = 0; i < kNumQueries; ++i) {
    pll.QueryDistance(vs[i], ws[i]);
  }
  double elapsed_time = GetCurrentTimeSec() - time_start;
  cout << "average query time: "
       << elapsed_time / kNumQueries * 1E6
       << " microseconds" << endl;

  exit(EXIT_SUCCESS);
}

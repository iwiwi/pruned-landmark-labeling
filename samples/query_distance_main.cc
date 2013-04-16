#include <cstdlib>
#include <iostream>
#include "pruned_landmark_labeling.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: construct_index INDEX" << endl;
    exit(EXIT_FAILURE);
  }

  PrunedLandmarkLabeling<> pll;
  if (!pll.LoadIndex(argv[1])) {
    cerr << "error: Load failed" << endl;
    exit(EXIT_FAILURE);
  }

  for (int u, v; cin >> u >> v; ) {
    cout << pll.QueryDistance(u, v) << endl;
  }
  exit(EXIT_SUCCESS);
}

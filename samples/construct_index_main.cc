#include <cstdlib>
#include <iostream>
#include "pruned_landmark_labeling.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "usage: construct_index GRAPH INDEX" << endl;
    exit(EXIT_FAILURE);
  }

  PrunedLandmarkLabeling<> pll;
  if (!pll.ConstructIndex(argv[1])) {
    cerr << "error: Load or construction failed" << endl;
    exit(EXIT_FAILURE);
  }
  pll.PrintStatistics();

  if (!pll.StoreIndex(argv[2])) {
    cerr << "error: Store failed" << endl;
    exit(EXIT_FAILURE);
  }
  exit(EXIT_SUCCESS);
}

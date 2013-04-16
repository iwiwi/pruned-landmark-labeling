#include <cstdio>
#include <cstdlib>
#include "pruned_landmark_labeling.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    fprintf(stderr, "usage: construct_index GRAPH INDEX\n");
    exit(EXIT_FAILURE);
  }

  PrunedLandmarkLabeling<> pll;
  if (!pll.ConstructIndex(argv[1])) {
    fprintf(stderr, "error: Load or construction failed\n");
    exit(EXIT_FAILURE);
  }
  if (!pll.StoreIndex(argv[2])) {
    fprintf(stderr, "error: Store failed\n");
    exit(EXIT_FAILURE);
  }
  exit(EXIT_SUCCESS);
}

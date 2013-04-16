Pruned Landmark Labeling
========================

Pruned landmark labeling is a new shortest-path distance querying index for real-world graphs, such as social networks, web graphs, biological networks and computer networks.

## Usage

### From CUI Interface

    $ make
    $ bin/construct_index samples/graph_example.tsv index_file
    $ bin/query_distance index_file <<< "1 4"
    2

* Execute `make` to build programs.
* Execute `bin/construct_index` to construct an index from a graph.
* Execute `bin/query_distance` and write pairs of vertices to STDIN to query distance between pairs of vertices.


### From Your Program

    PrunedLandmarkLabeling<> pll;
    pll.ConstructIndex(graph_file);
    cout << pll.QueryDistance(0, 1) << endl;

* Call `ConstructIndex` to construct an index from a graph.
* Call `QueryDistance` to query distance between two vertices.
* Call `LoadIndex` and `StoreIndex` to load and store an index.

For further information, please see `pruned_landmark_labeling.h`, samples and tests.

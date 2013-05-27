Pruned Landmark Labeling
========================

Pruned landmark labeling is a new shortest-path distance querying algorithm for real-world graphs, such as social networks, web graphs, biological networks and computer networks.

## Advantages
The algorithm has the following advantages (for details, please see our paper):

* **Fast** --- it answers distance queries in microseconds,
* **Scalable** --- it can be applied to networks with hundreds of millions of edges,
* **Exact** --- unlike approximate methods, it always answers exactly correct distance, and
* **Almost Parameter Free** --- unlike other state-of-the-art methods, it does not require any parameter tuning.

Moreover, this implementation is:

* **Easy to Use** --- by copying only one header file to your project, you can start using the index.

## Usage
Given a graph, it first constructs an index. Then, using the index, it can quickly answer distance between two vertices.

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
    pll.ConstructIndex(edge_list);
    cout << pll.QueryDistance(1, 4) << endl;

* Call `ConstructIndex` to construct an index from a graph (an edge list or a file).
* Call `QueryDistance` to query distance between two vertices.
* Call `LoadIndex` and `StoreIndex` to load and store an index.

For further information, please see `pruned_landmark_labeling.h`, samples and tests.

### Details

* In a graph file, each line should contain two integers describing an edge (see `samples/graph_example.tsv`).
* Vertices should be described by integers starting from zero.
* Program `bin/query_distance` reads pairs of vertices until EOF, thus you can use it to process multiple pairs of vertices at once.
* Execute `make test` to run tests (*google-gtest* is required).

## References

* Takuya Akiba, Yoichi Iwata, and Yuichi Yoshida, **[Fast Exact Shortest-Path Distance Queries on Large Networks by Pruned Landmark Labeling](http://www-imai.is.s.u-tokyo.ac.jp/~takiba/papers/sigmod13_pll.pdf)**.
In *SIGMOD 2013*, to appear.

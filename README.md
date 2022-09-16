# kPlexS: Maximum k-Plex Computation

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "kPlexS", which corresponds to the kPlexS algorithm.

For generating the algorithm kPlexF, please define NO_TRUSS_PRUNE in Utility.h by adding the following line and then recompiling the code.
```
#define NO_TRUSS_PRUNE
```

For generating the algorithm BnB-ct, please define NAIVE in Utility.h by adding the following line and then recompiling the code.
```
#define NAIVE
```

## Run the code

```sh
$ ./kPlexS -g {path_to_graph} -a exact -k {k_value} -o
```

An example of computing the exact maximum 3-plex for the dataset CA-GrQc is as follows
```sh
$ ./kPlexS -g datasets/CA-GrQc -a exact -k 2 -o
```

## Data format
Two data formats are supported. The default data format is "edges.txt", which contains a list of undirected edges represented as vertex pairs. The first line contains two numbers, representing the number of vertices and the number of undirected edges.

The more time-efficient format is the binary format; please replace "graph->read_graph();" by "graph->read_graph_binary();" in main.cpp. Each graph is represented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin). More details of the data format can be found in [https://lijunchang.github.io/Cohesive_subgraph_book/datasets](https://lijunchang.github.io/Cohesive_subgraph_book/datasets)

# kPlexS: Maximum k-Plex Computation

This repository implements the maximum k-plex computation algorithm proposed in our VLDB 2023 paper. If you are using the code, please cite our paper.
<pre>
Lijun Chang, Mouyi Xu and Darren Strash.
<a href="https://lijunchang.github.io/pdf/2022-Maximum-kPlex.pdf">Efficient Maximum k-Plex Computation over Large Sparse Graphs.</a>
Proc. VLDB Endow. 16(2), (2022)
</pre>

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "kPlexS", which corresponds to the kPlexS algorithm.

## Run the code

```sh
$ ./kPlexS -g {path_to_graph} -a exact -k {k_value} -o
```

An example of computing the exact maximum 3-plex for the dataset CA-GrQc is as follows
```sh
$ ./kPlexS -g datasets/CA-GrQc -a exact -k 2 -o
```

## Data format
Two data formats are supported. The default data format is "edges.txt", which contains a list of undirected edges represented as vertex pairs. The first line contains two numbers n and m, representing the number of vertices and the number of undirected edges, respectively. Note that, the vertex ids must be between 0 and n-1.

The more time-efficient format is the binary format; to read the input graph from this format, please add "-b" when running the code. Each graph is represented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin). More details of the data format can be found in [https://lijunchang.github.io/Cohesive_subgraph_book/datasets](https://lijunchang.github.io/Cohesive_subgraph_book/datasets)

## Algorithms kPlexF and BnB-ct

For generating the algorithm kPlexF, please undefine _SECOND_ORDER_PRUNING_ in KPlex_BB_matrix.h by commenting (or removing) the following line and then recompiling the code.
```
// #define _SECOND_ORDER_PRUNING_
```

For generating the algorithm BnB-ct, please define NAIVE in Utility.h by adding the following line and then recompiling the code.
```
#define NAIVE
```

# kPlexS: Maximum k-Plex Computation

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "kPlexS".

## Run the code

Different algorithms can be invoked by executing "MC-BRB".
```sh
$ ./MC-BRB {alg} {graph_directory}
```
alg is chosen from "MC-DD, MC-EGO, MC-BRB, verify"

An example of computing the exact maximum clique for the dataset CA-GrQc is as follows
```sh
$ ./MC-BRB MC-BRB datasets/CA-GrQc
```

## Data format
Each graph is represented by two binary files, b_adj.bin and b_degree.bin (e.g. datasets/CA-GrQc/b_adj.bin and datasets/CA-GrQc/b_degree.bin). More details of the data format can be found in [https://lijunchang.github.io/Cohesive_subgraph_book/datasets](https://lijunchang.github.io/Cohesive_subgraph_book/datasets)

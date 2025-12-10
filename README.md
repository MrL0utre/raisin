# Raisin project

## Short description

Raisin has been designed to solve the problem of circuit placement onto multi-FPGA architectures. It models the circuit to map as a set of red-black, directed, acyclic hypergraphs (DAHs). Hypergraph vertices can be either red vertices (which represent registers and external I/O ports) or black vertices (which represent internal combinatorial circuits). Vertices bear multiple weights, which define the types of resources needed to map the circuit (e.g., registers, ALUs, etc.). Every hyper-arc comprises a unique source vertex, all other ends of the hyper-arcs being sinks (which models the transmission of signals through circuit wiring). A circuit is consequently represented as set of DAHs that share some of their red vertices.

Target architectures are described by their number of target parts, the maximum resource capacities within each target part, and the connectivity between target parts.

The main metric to minimize is the length of the longest path between two red vertices, that is, the critical path that signals have to traverse
during a circuit compute cycle, which correlates to the maximum frequency at which the circuit can operate on the given target architecture.

Raisin computes a partition in which resource capacity constraints are respected and the critical path length is kept as small as possible, while reducing the number of cut hyper-arcs.
It produces an assignment list, which describes, for each vertex of the hypergraphs, the part to which the vertex is assigned.

## Features

- Partitioning 
- Clustering
- Refinement of a partition

## Installation

Run the command ```make``` to create the executable named *raisin* in the folder: *exe*. 
```sh
make 
```
To create executables for tests, you should run the command ```make``` followed by:
- **test_rbh**, for testing functions : load, write and convert red-black hypergraph.
- **test_crbh**, for testing clustering and coarsening algorithms.
- **test_ipart**, for testing initial partitioning algorithms: *Derived Breadth-First Search*; *Derived Depth-First Search*; *Critical Component Partitioning*.
- **test_fm**, for testing refinement algorithm taking target topology into account.
- **test_m**, for testing multilevel scheme (not stable).


To run clustering algorithm, you should select the *cluster* mode, and specify the algorithm: hem for heavy edge matching, and bsc for binary search clustering [2]. By default, the architecture of the target topology is a fully connected target topology.

```sh
raisin path_to_graph_file cluster algorithm(hem,bsc) size integer [bfactor integer] [partfile string] [archfile string]
```

To run partitioning algorithm, you should select the *part* mode, and specify the algorithm: dbfs for derived breadth-first search, ddfs for derived depth-first search and ccp for critical component partitioning [1,3].

```sh
raisin path_to_graph_file part algorithm(dbfs,ddfs,ccp) part_number integer [bfactor integer] [partfile string] [archfile string]
```

To run refinement algorithm, you should select the *refine* mode, and specify the algorithm: kfm for the classical K-way Fiduccia-Mattheyses refinement, and dkfm for Delay K-way Fiduccia-Mattheyses [1].

```sh
raisin path_to_graph_file refine algorithm(kfm,dkfm) part_number integer [bfactor integer] [perform integer] [tolerance integer] [partfile string] [archfile string]
```

## Development

Want to contribute? Great!


## References
[1] : Julien Rodriguez, François Galea, François Pellegrini et Lilia Zaourar. “A Hypergraph Model and Associated Optimization Strategies for Path Length-Driven Netlist Partitioning”. In : ICCS 2023 - 23rd International Conference on Computational Science. Sous la dir. de Jiří Mikyška, Clélia de Mulatier, Maciej Paszynski, Valeria V. Krzhizhanovskaya, Jack J. Dongarra et Peter M.A. Sloot. T. 10475. Lecture Notes in Computer Science. Prague, Czech Republic : Springer, juill. 2023, p. 652-660. doi : 10.1007/978-3-031-36024-4\_50. url : https://hal.science/hal-04379716.

[2] : Julien Rodriguez, François Galea, François Pellegrini et Lilia Zaourar. “Hypergraph Clustering with Path-Length Awareness”. In : ICCS 2024 - 24th International Conference on Computational Science. T. 14836. Lecture Notes in Computer Science. Malaga, Spain, juill. 2024, p. 90-104. doi : 10.1007/978-3-031-63775-9\_7. url : https://hal.science/hal-04706759.

[3] : Julien Rodriguez. Circuit partitioning for multi-FPGA platforms. Theses. Université de Bordeaux, sept. 2024. url : https://theses.hal.science/tel-04731886.


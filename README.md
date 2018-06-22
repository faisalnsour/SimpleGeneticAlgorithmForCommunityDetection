# A Simple Genetic Algorithm For Community Detection
Course project - a genetic algorithm that finds (near) optimal community structure for a network.

This GA uses locus-based adjacency (LBA) to represent a network and its (discrete) community structure. Consider a small graph, nodes labeled 0-4. One possible community structure can be captured with this LBA (i is index, n is node):

| i | n |
|---|---|
| 0 | 1 |
| 1 | 2 |
| 2 | 0 |
| 3 | 4 |
| 4 | 3 |

The "chain" of edges captured here represent a community. So nodes 0, 1, and 2 form a community, and nodes 3 and 4 form another. Using this representation allows the GA to its usual magic of selection, crossover, and mutation.

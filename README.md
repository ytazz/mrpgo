# mrpgo (Multi-Resolution Pose-Graph Optimizer)

## Abstract

A pose-graph optimization algorithm based on multi-resolution transform and parallel block-Gauss-Seidel iteration.
Check out the [RA-L paper](http://www2.kobe-u.ac.jp/~tazaki/docs/tazaki_ral2022.pdf) for details.

## Building

Use CMake.

### Dependent libraries

- Eigen3
- Intel MKL
- Suitesparse

## Running

mrpgo_app is a command-line application.
Run it without argument to display help.
It accepts input files in g2o format.
It recognizes VERTEX_SE2/EDGE_SE2 for 2D posegraphs,
 and VERTEX_SE3/EDGE_SE3:QUAT for 3D posegraphs.
Output file can be specified with -o option.

Visualization capability is not implemented.
Load the output posegraph file on g2o_viewer (for example) to see if its properly optimized.


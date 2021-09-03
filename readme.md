A collection of example GIS algorithms.

At the moment this includes:
- Shortest path (Djikstra)
- Eulerian / Tree tests for networks
- Bounding-box computation
- Delaunay triangulation (Sweephull)

There's also a partial implementation of Voronoi diagram generation using Fortune's algorithm, but it's a) only partially complete and b) known to be buggy.

A command-line program is included that will generate Graphviz and Gnuplot files for visualisation.

I've also included a shell script `diagrams.fish` for generating images from these, but unless you've got Fish installed you probably won't be able to use it. Should be easy enough to adapt to bash or similar though.
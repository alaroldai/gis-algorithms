A collection of example GIS algorithms.

At the moment this includes:
- Shortest path (Djikstra)
- Eulerian / Tree tests for networks
- Bounding-box computation


A command-line program is included that will generate Graphviz files for visualisation. You can produce pdfs from these with something like

```fish
for f in gen/*.dot; dot -Tpng $f > $f.png; end
```
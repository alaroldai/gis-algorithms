# On Delaunay Triangulation

The algorighm used here is an implementation of the "sweephull" or "S-Hull" algorithm by David Sinclair, a description of which can be found at [s-hull.org](https://www.s-hull.org/paper/s_hull.pdf).

In brief, the algorithm:
- Pick one of the input points as a "pivot" point and sort the points by
  distance to the pivot
- Construct an initial "seed" triangle from
  - the pivot,
  - it's nearest neighbour, and
  - the point that forms the smallest circumcircle with those two points
- Process the other points in order from closest to farthest away from the
  "pivot" point, keeping track of a convex hull of the set of points that have
  been added. This is the titular "sweephull".
  Each new point will form a triangle that shares one of the edges with the
  convex hull
- Track the triangles as they are created, including tracking neighbor
  relationships between triangles. These relationships form a graph that is
  used in the final step.
- Once all the points have been processed, the triangles form a valid
  triangulation but not a Delaunay triangulation.
  The last step is to process each pair of neighbouring triangles, check if
  they're a valid Delaunay pair (i.e., no point falls inside the circumcircle of
  the other three), and if not, "flip" the shared edge.
  After this is finished, the triangulation is a Delaunay triangulation.

Complexity: O(n*log(n)) is claimed by the S-Hull paper cited above.
As best I can tell, this is limited by the sort-by-radial distance at the
beginning of the algorithm, and the need to process the convex hull every time a
new point is added.

Implementation is in three parts:
- Initialisation (`new`) does the sort and constructs the seed triangle
- A `step` function implements the processing of a single point from the
  sorted point array. Iteration variables are stored in this struct.
- After the final point has been processed, `flip_step` does the final
  flipping pass over all the triangles.
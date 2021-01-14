To Do List
==========

Done
====

- Develop the fundamental domain, i.e., assign an :class:`IdealPoint` to each vertex of each simplex.

  - this is implemented in ``IdealPointFundamentalDomainVertexEngine``, briefly tested and the vertices are the same than what SnapPea kernel produces.

- Compute matrix taking three ideal points to three ideal points.

  - this is implemented in ``matrix_taking_triple_to_triple``.
  - by solving system of the three linear equations
  - thus different from how SnapPy does it in ``compute_matrices in snap/generators.py

- Implement an interval tree which we can use to efficiently look up tiles

- Compute center of an ideal tetrahedron and radius of maximal ball

  - Write helper which does this for the tetrahedron with vertices, 0, inf, 1, z.

Todo
====

- Double check the code in ``radius_and_center_of_largest_ball``

- Compute sinew

  - Use ``RealCuspCrossSection`` to compute a cross section
  - use the cross section to determine a center for each edge of the triangulation
  - compute center for each face of the triangulation
  - compute radius of sinew

- Compute distance between a finite and an ideal triangle

  - compute distance between the two planes containing the triangles and the end points of the shortest arcs between the planes, check whether the end points are in the triangles
  - compute distance between line and plane and the end points of the shortest arc and cheeck whether in the line segment and in the triangle
  - similar line and point
  - similar line and line (using Maria's trick)

- Using above, given a lambda, decide whether a face of a translate of a fundamental domain has distance at least lambda from sinew, i.e., whether we do not need to tile further

  - Implement early bail-outs. For example, if the distane between the given face and the center of the sinew is less than lambda, we can return False. If the distance between the given face and the center of the sinew is greater than lambda plus the radius of the sinew, we can return True.

- Implement the tiling algorithm

  - use data structure to lookup a tile

- Filter out powers and conjugates

- Ambitious: drill out geodesics



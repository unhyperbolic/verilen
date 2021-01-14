from snappy.snap.t3mlite import simplex

from snappy.verify.mathHelpers import interval_aware_min

from tilingEngineBase import TilingEngineBase

from upperHalfspace.projectivePoint import *
from upperHalfspace.idealTriangle import *
from upperHalfspace.finiteIdealTriangleDistanceEngine import *

from spineEngine import *

from sage.all import matrix

import heapq

__all__ = ['SpineTilingEngine']

class SpineTilingEngine(TilingEngineBase):
    """
    
    A tiling algorithm that is using a geometric spine to determine
    how far to tile out.

    More precisely, given a cut-off length lambda, it tiles until
    the distance of any unglued face to the geometric spine is larger
    than lambda.

    Note: This is not finished yet!!!
    For now, it just keeps tiling every time you call add_next_tile
    which just adds a random tile.

    >>> from testTilingEngine import *
    >>> from snappy import *
    >>> M = Manifold("s776")
    >>> M = Manifold("o9_24124")
    >>> shapes = M.verify_hyperbolicity(bits_prec = 512)[1]
    >>> t = SpineTilingEngine.fromManifoldAndShapesMatchingSnapPea(M, shapes)

    >>> tile = t.add_tile(t.mcomplex.GeneratorMatrices[1])

    >>> bool(t.find_tile(t.mcomplex.GeneratorMatrices[1]))
    True
    >>> bool(t.find_tile(t.mcomplex.GeneratorMatrices[2]))
    False

    >>> for i in range(300):
    ...     tile = t.add_next_tile()
    
    >>> simple_tiling_engine_consistency_check(t)

    """


    @staticmethod
    def fromManifoldAndShapesMatchingSnapPea(snappyManifold, shapes):
        e = SpineEngine.fromManifoldAndShapesMatchingSnapPea(
            snappyManifold, shapes)
        return SpineTilingEngine(e.mcomplex)

    def __init__(self, mcomplex):
        # A list of pairs (lower bound, tile, generator) of corresponding to
        # unglued faces.

        self.lower_bounds_and_unglued_generators = []

        super(SpineTilingEngine, self).__init__(mcomplex)

    def compute_distance_tile_face_to_spine(self, tile, corner):
        tet = self.mcomplex.Tetrahedra[0]
        CIF = tet.ShapeParameters[simplex.E01].parent()

        vertices = [
            ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(
                CIF,
                tile.vertexToIdealPoint[corner.Tetrahedron.Class[V]])
            for V in simplex.ZeroSubsimplices if V & corner.Subsimplex ]

        ideal_triangle = IdealTriangle(vertices)

        return interval_aware_min(
            FiniteIdealTriangleDistanceEngine(
                finite_triangle, ideal_triangle).compute_lowerbound()
            for tet in self.mcomplex
            for finite_triangle in tet.spineTriangles.values())

    def compute_distance_tile_gen_to_spine(self, tile, g):
        return interval_aware_min(
            self.compute_distance_tile_face_to_spine(
                tile, corner)
            for (corner, other_corner), perm in self.mcomplex.Generators[g])

    def add_initial_tile(self):
        tet = self.mcomplex.Tetrahedra[0]
        CIF = tet.ShapeParameters[simplex.E01].parent()
        m = matrix.identity(CIF, 2)
        self.add_tile(m, is_initial = True)

    def add_tile(self, m, is_initial = False):
        """
        Create and add tile corresponding to matrix m.

        Keep track of the unglued faces.
        """

        tile = self.create_tile_and_add_to_tree(m)
        
        unglued_generators = self.process_tile(tile)

        for g in unglued_generators:
            entry = (self.compute_distance_tile_gen_to_spine(tile, g).lower(),
                     tile,
                     g)

            heapq.heappush(self.lower_bounds_and_unglued_generators, entry)

        print(self.lower_bounds_and_unglued_generators)

        return tile

    def add_next_tile(self):
        """
        Find the next unglued face and create the tile neighboring
        along that unglued face.
        """

        while True:
            lower_bound, tile, g = heapq.heappop(self.lower_bounds_and_unglued_generators)
            if not tile.generatorToNeighboringTile.has_key(g):
                # TODO!
                # Add an additional if here that skips if all the faces
                # corresponding to (tile, generator) have distane
                # greater than the cut-off length lambda from the geometric
                # spine in the fundamental polyhedron.
                #
                # This can be done using the FiniteIdealTriangleDistanceEngine.
                #
                m = tile.matrix * self.mcomplex.GeneratorMatrices[g]
                return self.add_tile(m)

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()

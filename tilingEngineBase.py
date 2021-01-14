from snappy.snap.mcomplex_base import *

from snappy.snap.t3mlite import simplex

from snappy.verify.interval_tree import IntervalTree

from snappy.verify.upper_halfspace.ideal_point import *

from sage.all import matrix

class TilingEngineBase(McomplexEngine):
    """
    Base class for algorithms tiling H^3 with fundamental polyhedra F.
    Here, a tile is some data structure that stores the matrix m such
    that the tile is given by mF. A tile also stores pointers to the
    neighboring tiles.
    
    The class needs to be instantiated with an Mcomplex object that is a
    fundamental polyhedron.

    The basic function of this class is to keep track of what tiles have
    already been enumerated and update the adjacency relations between the
    tiles when adding a tile.

    To do this, the class will choose a base point p for the fundamental
    polyhedron.
    """


    class Tile(object):
        def __init__(self, matrix, center):
            """
            A tile which is the image of the fundamental polyhedron
            under the given matrix (a SageMath matrix with complex interval
            field entries). center needs to be an instance of ``FinitePoint``
            and needs to be image of the base point of the fundamental
            polyhedron by the given matrix.
            """

            self.matrix = matrix
            self.center = center

            # Give the point (type IdealPoint) in CP^1 (regarded as
            # one-point compactification) where a vertex of this
            # tile is located at.
            # The key is the corresponding vertex in the fundamental
            # polyhedron (of type snappy.snap.t3mlite.Vertex).
            self.vertexToIdealPoint = {}

            # Pointers to the neighboring tiles.
            # The key is the generator of the face-pairing representation
            # of the fundamental group that the SnapPea kernel is using
            # (an integer with the usual convention that -i is the inverse of
            # i).
            self.generatorToNeighboringTile = {}

    def __init__(self, mcomplex):
        """
        Construct from an Mcomplex that was created with
        FundamentalPolyhedronEngine.from_manifold_and_shapes_matching_snappea.
        """

        super(TilingEngineBase, self).__init__(mcomplex)

        # An interval tree to look-up existing tiles by translated base point.
        self.intervalTree = IntervalTree()

        # Compute a base point of the fundamental polyhedron
        tet = self.mcomplex.ChooseGenInitialTet
        self.baseTetInRadius, self.baseTetInCenter = compute_inradius_and_incenter(
            [ tet.Class[v].IdealPoint for v in simplex.ZeroSubsimplices])

        # Add a tile corresponding to the identity.
        self._add_initial_tile()

    def are_same_tile(self, center1, center2):
        """
        Let the center of a tile be the image of the base point in the
        fundamental polyhedron by the transformation taking the tile
        to the fundamental polyhedron.

        Given two objects of type ``FinitePoint`` that are intervals for
        such centers, decide whether they are for the same tile or not.

        This will either return True or False or fail if the intervals
        are too large to verify the answer.
        
        """

        dist = center1.dist(center2)

        # I think this can be 2 * R
        if dist < self.baseTetInRadius:
            return True
        # And this can be 0.
        if dist > self.baseTetInRadius:
            return False

        raise Exception(
            "Distance between two given centers of tiles cannot be verified "
            "to be small enough to be the same or large enough to be two different "
            "tiles")
        
    def find_tile(self, m):
        """
        Given a matrix m taking the fundamental polyhedron to a tile
        (a SageMath matrix with ComplexIntervalField entries), determine
        whether the tile has already been enumerated earlier.

        Return that tile or return None.
        """
        
        center = self.baseTetInCenter.translate_PGL(m)
        for tile in self.intervalTree.find(center.key_interval()):
            if self.are_same_tile(center, tile.center):
                return tile

        return None

    def create_tile_and_add_to_tree(self, m):
        """
        Given a matrix m taking the fundamental polyhedron to a tile,
        create the tile and add it to the acceleration structure for
        later look-up with find_tile.

        This returns the new tile. It does not update adjacency relations.
        """
        
        center = self.baseTetInCenter.translate_PGL(m)
        tile = TilingEngineBase.Tile(m, center)
        self.intervalTree.insert(center.key_interval(), tile)

        return tile

    def process_tile(self, tile):
        """
        Given a new tile (created with create_tile_and_add_to_tree), update
        adjacency relations for all existing tiles and compute the
        vertices of the new tile.

        Returns the list of unglued generators. This is a list of
        generators referring to the face-pairing presentation of the
        fundamental group (as integers with the negative referring to
        the inverse). A generator is in this list, if it takes the new
        tile to a tile that has not been enumerated previously yet.

        In other words, these are all the generators that need to be applied
        to the new tile to further the tiling.
        """

        unglued_generators = []

        for g, gen_m in sorted(self.mcomplex.GeneratorMatrices.items()):
            other_m = tile.matrix * gen_m
            
            other_tile = self.find_tile(other_m)
            if other_tile:
                tile.generatorToNeighboringTile[g] = other_tile
                other_tile.generatorToNeighboringTile[-g] = tile

                for (corner, other_corner), perm in self.mcomplex.Generators[g]:
                    for v in simplex.VerticesOfFaceCounterclockwise[corner.Subsimplex]:
                        vertex       = corner.Tetrahedron.Class[v]
                        other_vertex = other_corner.Tetrahedron.Class[perm.image(v)]
                        
                        tile.vertexToIdealPoint[vertex] = (
                                other_tile.vertexToIdealPoint[other_vertex])    
            else:
                unglued_generators.append(g)
                
        for vertex in self.mcomplex.Vertices:
            if not vertex in tile.vertexToIdealPoint:
                tile.vertexToIdealPoint[vertex] = apply_Moebius(
                    tile.matrix, vertex.IdealPoint)
    
        return unglued_generators

    def _add_initial_tile(self):
        tet = self.mcomplex.Tetrahedra[0]
        CIF = tet.ShapeParameters[simplex.E01].parent()
        m = matrix.identity(CIF, 2)
        return self.add_tile(m)

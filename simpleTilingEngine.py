from snappy.snap.fundamental_polyhedron import *

from upperHalfspace import *
from snappy.verify.upper_halfspace.ideal_point import *

from tilingEngineBase import TilingEngineBase

from collections import deque

class SimpleTilingEngine(TilingEngineBase):
    """
    
    >>> from testTilingEngine import *
    >>> from snappy import *
    >>> M = Manifold("s776")
    >>> t = SimpleTilingEngine.from_manifold_and_shapes(M, M.verify_hyperbolicity()[1])

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
    def from_manifold_and_shapes(
            manifold, shapes, normalize_matrices = False, match_kernel = True):
        f = FundamentalPolyhedronEngine.from_manifold_and_shapes(
            manifold, shapes,
            normalize_matrices = normalize_matrices,
            match_kernel = match_kernel)
        m = f.mcomplex
        t = SimpleTilingEngine(m)

        return t

    def __init__(self, mcomplex):
        self.unglued_generators = deque()

        super(SimpleTilingEngine, self).__init__(mcomplex)

    def add_tile(self, m):
        tile = self.create_tile_and_add_to_tree(m)
        
        unglued_generators = self.process_tile(tile)

        for g in unglued_generators:
            self.unglued_generators.append((tile, g))

        return tile

    def add_next_tile(self):
        while True:
            tile, g = self.unglued_generators.popleft()
            if not g in tile.generatorToNeighboringTile:
                m = tile.matrix * self.mcomplex.GeneratorMatrices[g]
                return self.add_tile(m)
                
def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()

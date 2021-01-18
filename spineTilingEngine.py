from snappy.snap.t3mlite import simplex

from snappy.verify.mathHelpers import interval_aware_min
from snappy.verify.upper_halfspace.ideal_point import _adjoint2
from snappy.verify.upper_halfspace.finite_point import *

from tilingEngineBase import TilingEngineBase

from upperHalfspace.projectivePoint import *

from spineEngine import *

from sage.all import matrix, sqrt, RealIntervalField, Infinity, acosh

# __all__ = ['SpineTilingEngine']

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
    >>> t = SpineTilingEngine.from_manifold_and_shapes(M, shapes)
    >>> t.target_radius = shapes[0].real().parent()(0.1)

    >>> tile = t.add_tile(t.mcomplex.GeneratorMatrices[1])

    >>> bool(t.find_tile(t.mcomplex.GeneratorMatrices[1]))
    True
    >>> bool(t.find_tile(t.mcomplex.GeneratorMatrices[2]))
    False

    >>> for i in range(100):
    ...     tile = t.add_next_tile()
    ...     print(i, bool(tile))
    
    >>> simple_tiling_engine_consistency_check(t)

    """


    @staticmethod
    def from_manifold_and_shapes(
            manifold, shapes, normalize_matrices = False, match_kernel = True):
        e = SpineEngine.from_manifold_and_shapes(
            manifold, shapes,
            normalize_matrices = normalize_matrices,
            match_kernel = match_kernel)
        return SpineTilingEngine(e.mcomplex)

    def __init__(self, mcomplex):
        # A list of pairs (lower bound, tile, generator) of corresponding to
        # unglued faces.

        self.unglued_generators = []

        super(SpineTilingEngine, self).__init__(mcomplex)

        z = self.mcomplex.Tetrahedra[0].ShapeParameters[simplex.E01]
        self.CIF = z.parent()
        self.RIF = z.real().parent()

        for tet in self.mcomplex.Tetrahedra:
            tet.SpineRadius = self.RIF(0)

            for f in simplex.TwoSubsimplices:
                tet.SpineRadius = tet.SpineRadius.max(
                    tet.InCenter.dist(tet.Class[f].InCenter))
                
            for e in simplex.OneSubsimplices:
                tet.SpineRadius = tet.SpineRadius.max(
                    tet.InCenter.dist(tet.Class[e].MidPoint))
            
    def finish_tiling(self):
        while self.add_next_tile():
            pass

    def needs_to_be_added(self, m):
        
        for tet in self.mcomplex.Tetrahedra:
            for f in simplex.TwoSubsimplices:
                if not tet.Neighbor[f]:

                    projectivePoints = [
                        ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(
                            self.CIF, tet.Class[v].IdealPoint).translate(m)
                        for v in simplex.ZeroSubsimplices
                        if v & f ]
                    
                    t = _adjoint2(
                        ProjectivePoint.matrix_taking_0_1_inf_to_given_points(
                            *projectivePoints))

                    for tet1 in self.mcomplex.Tetrahedra:
                        if not has_distance_larger(
                                tet1.InCenter.translate_PGL(t),
                                tet1.SpineRadius + self.target_radius):
                            return True

        return False

    def add_tile(self, m, is_initial = False):
        """
        Create and add tile corresponding to matrix m.

        Keep track of the unglued faces.
        """

        tile, unglued_generators = self.create_tile(m)
        for g in unglued_generators:
            self.unglued_generators.append((tile, g))

        return tile

    def add_next_tile(self):
        """
        Find the next unglued face and create the tile neighboring
        along that unglued face.
        """

        while self.unglued_generators:
            tile, g = self.unglued_generators.pop()
            if not g in tile.generatorToNeighboringTile:
                m = tile.matrix * self.mcomplex.GeneratorMatrices[g]
                if self.needs_to_be_added(m):
                    return self.add_tile(m)

def has_distance_larger(finPoint, dist):
    CIF = finPoint.z.parent()

    t = sqrt(  finPoint.z.imag() ** 2
             + finPoint.t ** 2)
    pt = FinitePoint(CIF(finPoint.z.real(), 0), t)

    d = finPoint.dist(pt)

    if d > dist:
        return True

    m = matrix(CIF, [[0,-1],[1,-1]])

    debug = finPoint

    for i in range(3):
        if finPoint.z.real() < 0:
            t = sqrt(  finPoint.z.real() ** 2
                     + finPoint.z.imag() ** 2
                     + finPoint.t ** 2)
            pt = FinitePoint(CIF(0), t)
            d = finPoint.dist(pt)
            if d > dist:
                return True

        finPoint = finPoint.translate_PSL(m)

    return False

def length_spectrum_with_multiples(M, cut_off, bits_prec = 53):
    
    is_hyp, shapes = M.verify_hyperbolicity(bits_prec = bits_prec)
    if not is_hyp:
        raise ValueError("Manifold does not seem to be hyperbolic.")
    
    RIF = RealIntervalField(bits_prec)

    t = SpineTilingEngine.from_manifold_and_shapes(M, shapes)
    t.target_radius = RIF(cut_off)
    t.finish_tiling()

    tiles = t.intervalTree.find(RIF(-Infinity, Infinity))

    lengths = [ 2 * acosh(tile.matrix.trace() / sqrt(tile.matrix.det()) / 2) for tile in tiles ]

    return [ length for length in lengths if length.real() < t.target_radius ]

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()

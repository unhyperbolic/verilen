from snappy.snap.t3mlite import simplex

from snappy.verify.mathHelpers import interval_aware_min
from snappy.verify.upper_halfspace.ideal_point import _adjoint2
from snappy.verify.upper_halfspace.finite_point import FinitePoint, _abs_sqr

from tilingEngineBase import TilingEngineBase

from upperHalfspace.projectivePoint import *
from upperHalfspace.fixedPoints import *

from spineEngine import *

from sage.all import matrix, sqrt, exp, RealIntervalField, Infinity, acosh

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
            manifold, shapes, cut_off, normalize_matrices = True, match_kernel = True):
        e = SpineEngine.from_manifold_and_shapes(
            manifold, shapes,
            normalize_matrices = normalize_matrices,
            match_kernel = match_kernel)
        return SpineTilingEngine(e.mcomplex, cut_off)

    def __init__(self, mcomplex, cut_off):
        # A list of pairs (lower bound, tile, generator) of corresponding to
        # unglued faces.

        self.unglued_generators = []

        super(SpineTilingEngine, self).__init__(mcomplex)

        z = self.mcomplex.Tetrahedra[0].ShapeParameters[simplex.E01]
        self.CIF = z.parent()
        self.RIF = z.real().parent()
        self.cut_off = cut_off

        for tet in self.mcomplex.Tetrahedra:
            tet.SpineRadius = self.RIF(0)

            for f in simplex.TwoSubsimplices:
                tet.SpineRadius = tet.SpineRadius.max(
                    tet.InCenter.dist(tet.Class[f].InCenter))
                
            for e in simplex.OneSubsimplices:
                tet.SpineRadius = tet.SpineRadius.max(
                    tet.InCenter.dist(tet.Class[e].MidPoint))

            tet.SpineRadiusPlusCutOff = tet.SpineRadius + self.cut_off
            e = exp(tet.SpineRadiusPlusCutOff)
            tet.CutOffFactor = (e * e - 1) / (2 * e)

        for tet in self.mcomplex.Tetrahedra:
            tet.transform_taking_face_to_0_1_inf = {
                f : _compute_transform_taking_face_to_0_1_inf(
                    tet, f, self.CIF)
                for f in simplex.TwoSubsimplices }
                
        self.initial_tile.word = []
            
    def finish_tiling(self):
        while self.add_next_tile():
            pass

    def needs_to_be_added(self, m):
        
        minv = _adjoint2(m)
        
        for tet in self.mcomplex.Tetrahedra:
            for f in simplex.TwoSubsimplices:
                if not tet.Neighbor[f]:

                    mt = tet.transform_taking_face_to_0_1_inf[f] * minv

                    for tet1 in self.mcomplex.Tetrahedra:
                        if not has_distance_larger(
                                tet1.InCenter.translate_PSL(mt),
                                tet1.CutOffFactor):
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
            tile, g = self.unglued_generators.pop(0)
            if not g in tile.generatorToNeighboringTile:
                m = tile.matrix * self.mcomplex.GeneratorMatrices[g]
                if self.needs_to_be_added(m):
                    new_tile = self.add_tile(m)
                    new_tile.word = tile.word + [ g ]
                    return new_tile

    def all_tiles(self):
        return self.intervalTree.find(self.RIF(-Infinity, Infinity))

def _compute_transform_taking_face_to_0_1_inf(tet, f, CIF):
    projectivePoints = [
        ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(
            CIF, tet.Class[v].IdealPoint)
        for v in simplex.ZeroSubsimplices
        if v & f ]

    m = ProjectivePoint.matrix_taking_0_1_inf_to_given_points(*projectivePoints)
    
    return _adjoint2(m / sqrt(m.det()))

def has_distance_larger(finPoint, distFactor):
    d = finPoint.t * distFactor

    if abs(finPoint.z.imag()) > d:
        return True

    if finPoint.z.real() < 0:
        return abs(finPoint.z) > d
    if finPoint.z.real() > 1:
        return abs(finPoint.z - 1) > d

    RIF = finPoint.z.real().parent()
    if (finPoint.z.real() - RIF(0.5)) ** 2 + finPoint.z.imag() ** 2 + finPoint.t ** 2 < RIF(0.25):
        f = _rotate_helper(finPoint)
        return abs(f.z) > f.t * distFactor

    return False

def _rotate_helper(finPoint):
    tSqr = finPoint.t ** 2
    num   = (finPoint.z - 1) * finPoint.z.conjugate() + tSqr
    denom = _abs_sqr(finPoint.z) + tSqr
    
    return FinitePoint(num / denom, finPoint.t / denom)

def is_not_parabolic(m, spineEngine):
    tr = m.trace()

    if tr != 2 and tr != -2:
        return True

    fixed_pts, bounding_halfsphere = (
        fixed_projective_points_or_bounding_halfsphere_for_matrix(m))

    if not bounding_halfsphere:
        raise Exception("Could not find bound for fixed point of presumambly parabolic matrix")
    
    for tet in spineEngine.mcomplex.Tetrahedra:
        if not bounding_halfsphere.is_outside(tet.InCenter):
            raise Exception("Could not bound away fixed points of presumambly parabolic matrix from tet incenter")

    for face in spineEngine.mcomplex.Faces:
        if not bounding_halfsphere.is_outside(face.InCenter):
            raise Exception("Could not bound away fixed points of presumambly parabolic matrix from face incenter")
        
    for edge in spineEngine.mcomplex.Edges:
        if not bounding_halfsphere.is_outside(edge.MidPoint):
            raise Exception("Could not bound away fixed points of presumambly parabolic matrix from edge midpoint")

    return False

def length_from_matrix(m):
     t = m.trace()
     if t.real().center() < 0:
         # Trace only determined up to sign
         # We prefer the one with positive real part because
         # SageMath's acosh has undesirable branching behavior.
         t = -t

     return 2 * acosh(t / 2)

def get_tiling_engine(M, cut_off, bits_prec = 53):
    is_hyp, shapes = M.verify_hyperbolicity(bits_prec = bits_prec)
    if not is_hyp:
        raise ValueError("Manifold does not seem to be hyperbolic.")
    
    RIF = RealIntervalField(bits_prec)

    t = SpineTilingEngine.from_manifold_and_shapes(M, shapes, RIF(cut_off))
    t.finish_tiling()

    return t

def length_spectrum_with_multiples(M, cut_off, bits_prec = 53):
    tiling_engine = get_tiling_engine(M, cut_off, bits_prec)

    lengths = [ length_from_matrix(tile.matrix)
                for tile in tiling_engine.all_tiles()
                if ((not tile is tiling_engine.initial_tile) and
                    is_not_parabolic(tile.matrix, tiling_engine)) ]

    return [ length for length in lengths
             if not length.real() > tiling_engine.cut_off ]

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()

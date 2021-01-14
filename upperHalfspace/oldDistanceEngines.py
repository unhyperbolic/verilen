from .projectivePoint import *
from .finitePoint import *
from .extendedMatrix import *
from .fixedPoints import *
from .reflections import *

from sage.all import sqrt, arccosh, arcsinh, log
import sage.all

class DistanceEngineBase(object):
    def __init__(self, *matrices):
        for dim, m in zip(self.dimensions, matrices):
            expected = - ((-1) ** dim)
            if ExtendedMatrix.get_orientation_sign(m) != expected:
                raise Exception(
                    "Expected reflection about %d subspace to act "
                    "by %d on orientation" % (dim, expected))

        self.matrices = matrices
        m = matrices[0] * matrices[1]
        if sum(self.dimensions) % 2 == 1:
            self.orientationPreservingComposition = m * m
        else:
            self.orientationPreservingComposition = m

        self.half_trace = ExtendedMatrix.get_trace_in_PSL(
            self.orientationPreservingComposition) / 2

    def fixed_projective_points_or_bounding_halfsphere(self):
        return fixed_projective_points_or_bounding_halfsphere_for_matrix(
            self.orientationPreservingComposition)

class PlanePlaneDistanceEngine(DistanceEngineBase):
    dimensions = [2, 2]

    def abs_half_real_trace(self):
        return abs(self.half_trace.real())

    def is_intersecting(self):
        return abs(self.half_real_trace()) < 1

    def is_not_intersecting_or_touching(self):
        return abs(self.half_real_trace()) > 1
    
    def is_not_identical(self):
        return _is_not_identity(self.orientationPreservingComposition)

    def dist(self):
        return arccosh(self.abs_half_real_trace())

    def endpoints(self):
        pts, max_height = self.fixed_projective_points_or_bounding_halfsphere()
        line = LineReflection.from_two_projective_points(*pts)
        
        return [ PointReflection.to_finite_point(m * line)
                 for m in self.matrices ]

class PlaneLineDistanceEngine(DistanceEngineBase):
    dimensions = [2, 1]

    def abs_half_real_trace(self):
        return abs(self.half_trace.real())

    def dist(self):
        return arccosh(self.abs_half_real_trace()) / 2        

    # More thinking required

class PlanePointDistanceEngine(DistanceEngineBase):
    dimensions = [2, 0]

    def half_imag_trace(self):
        return self.half_trace.imag()

    def dist(self):
        return arcsinh(self.half_imag_trace())

    def endpoint(self):
        pts, max_height = self.fixed_projective_points_or_bounding_halfsphere()
        line = LineReflection.from_two_projective_points(*pts)
        
        return PointReflection.to_finite_point(self.matrices[0] * line)

class LineLineDistanceEngine(DistanceEngineBase):
    dimension = [1, 1]

    def trace(self):
        return _trace(self.composition)

    def endpoint(self):
        return _intersection_of_two_lines_from_endpoints_and_matrix(
            self.endpoints_matrix1, self.composition)

class LinePointDistanceEngine(DistanceEngineBase):
    dimension = [1, 0]

################################################################################
#
# Various helper functions

def _intersection_of_two_lines_from_endpoints_and_matrix(pts_line1, matrix2):
    pts_line2, h = fixed_projective_points_or_bounding_sphere_for_matrix(matrix2)
    return _intersection_of_two_lines_from_endpoints(pts_line1, pts_line2)

def _intersection_of_two_lines_from_endpoints(pts_line1, pts_line2):
    m = matrix_taking_0_1_inf_to_given_points(
        pts_line1[0], pts_line1[1], pts_line2[0])

    z = pts_line2[1].translate(_adjoint2(m)).as_ideal_point()
    if z == Infinity:
        raise Exception("Unexpected Infinity")
    CIF = z.parent()
    zReal = z.real()
    t = (zReal * 2 - 1).sqrt() / 2
    
    return FinitePoint(CIF(zReal), t).translate_PGL(m)

def _is_not_identity(m):
    if isinstance(m, ExtendedMatrix):
        m = m.matrix
    return m[1,0] != 0 or m[0,1] != 0 or m[0,0] != m[1,1]

################################################################################
#
# TESTING

class _DistanceEngineTester(object):
    
    """
    A test rig for the derived classes of DistanceEngineBase.

    Run the test rig::

        sage: _DistanceEngineTester().run_tests()

    """

    def test_plane_plane(self):
        from sage.all import RealIntervalField, ComplexIntervalField
        CIF = ComplexIntervalField(120)
        RIF = RealIntervalField(120)

        pts = [
            [ ProjectivePoint(CIF(1)),
              ProjectivePoint(CIF(3)),
              ProjectivePoint(CIF(2,1)) ],
            [ ProjectivePoint(CIF(-1)),
              ProjectivePoint(CIF(-3)),
              ProjectivePoint(CIF(-2,1)) ] ]

        planes = [ 
            PlaneReflection.from_three_projective_points(*pt)
            for pt in pts ]

        p = PlanePlaneDistanceEngine(*planes)
        e = p.endpoints()
        if not abs(p.dist() - e[0].dist(e[1])) < 1e-20:
            raise Exception("Distance incorrect %r %r" % (
                    p.dist(), e[0].dist(e[1])))

        e_baseline = [
            FinitePoint(CIF( 1.5), RIF(0.75).sqrt()),
            FinitePoint(CIF(-1.5), RIF(0.75).sqrt()) ]

        for e0, e_baseline0 in zip(e, e_baseline):
            if not e0.dist(e_baseline0) < 1e-15:
                raise Exception("%s %s" % (e0, e_baseline0))
        
    def test_plane_line(self):
        from sage.all import CIF
        
        pts_plane = [ ProjectivePoint(CIF(-1)),
                      ProjectivePoint(CIF(1)),
                      ProjectivePoint(CIF(0,1)) ]
        pts_line = [ ProjectivePoint(CIF(sage.all.e)),
                     ProjectivePoint(CIF(-sage.all.e)) ]
        plane = PlaneReflection.from_three_projective_points(*pts_plane)
        line = LineReflection.from_two_projective_points(*pts_line)

        pp = PlaneLineDistanceEngine(plane, line)

        if not abs(1 - pp.dist()) < 1e-10:
            raise Exception("Wrong distance: %r" % pp.dist())
        

    def test_plane_point(self):
        from sage.all import RealIntervalField, ComplexIntervalField
        CIF = ComplexIntervalField(120)
        RIF = RealIntervalField(120)
        
        pts = [ ProjectivePoint(CIF(1.01)),
                ProjectivePoint(CIF(3)),
                ProjectivePoint(CIF(2,1)) ]

        pt = FinitePoint(CIF(2,5), RIF(3))

        plane = PlaneReflection.from_three_projective_points(*pts)
        
        p = PlanePointDistanceEngine(
            plane,
            PointReflection.from_finite_point(pt))
        e = p.endpoint()
        
        if not abs(p.dist() - pt.dist(e)) < 1e-20:
            raise Exception("Distance incorrect %r %r" % (
                    p.dist(), pt.dist(e)))

    def run_tests(self):
        
        self.test_plane_plane()
        self.test_plane_line()
        self.test_plane_point()

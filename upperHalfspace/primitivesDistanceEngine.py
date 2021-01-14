from snappy.verify.upper_halfspace.extended_matrix import *
from .fixedPoints import *

from sage.all import sqrt, arccosh, arcsinh, log
import sage.all

__all__ = ['PrimitivesDistanceEngine']

class PrimitivesDistanceEngine(object):
    """
    An engine to compute the distance between two primitives, i.e., plane,
    line, or point. Each primitive is represented by the involution whose
    fixed points fom that primitive.

    The constructor takes as first argument a list with the dimensions of the
    primitives and as second argument the involutions as
    :class:`ExtendedMatrix`\ s or SageMath ``matrix``::
    
         sage: from sage.all import *
         sage: planeReflection = ExtendedMatrix(matrix.identity(CIF, 2), True)
         sage: lineReflection = ExtendedMatrix(matrix(CIF, [[-3j, -4], [-2, 3j]]))
         sage: e = PrimitivesDistanceEngine([2, 1], [planeReflection, lineReflection])
         sage: e.dist()
         1.762747174039086?

    It works by computing the product of the two involutions, squaring it if
    necessary to create an orientation-preserving hyperbolic motion::
    
         sage: e.orientationPreservingComposition
         ExtendedMatrix([   17 -24*I]
         [ 12*I    17])

    And computing the trace of it (when normalized to SL(2,C)) and dividing it
    by two::

         sage: e.half_trace
         17
    """

    def __init__(self, dimensions, matrices):
        PrimitivesDistanceEngine._do_checks(dimensions, matrices)

        self.dimensions = dimensions
        m = matrices[0] * matrices[1]

        if sum(self.dimensions) % 2 == 1:
            self.orientationPreservingComposition = m * m
        else:
            self.orientationPreservingComposition = m

        self.half_trace = ExtendedMatrix.get_trace_in_PSL(
            self.orientationPreservingComposition) / 2

    @staticmethod
    def _do_checks(dimensions, matrices):
        if len(dimensions) != 2:
            raise Exception("Expected two dimensions")
        if len(matrices) != 2:
            raise Exception("Expected two matrices")

        for dim, m in zip(dimensions, matrices):
            expected = - ((-1) ** dim)
            if ExtendedMatrix.get_orientation_sign(m) != expected:
                raise Exception(
                    "Expected reflection about %d subspace to act "
                    "by %d on orientation" % (dim, expected))

    def is_stably_intersecting(self):
        if sum(self.dimensions) < 3:
            # Return false for all other cases
            return False

        # The Plane-plane and Plane-line case are distinguished in that they
        # can intersect such that they stay intersected under a small
        # perturbation to either.
        #
        # The motivation to detect this case is this: the geodesic fixed by the
        # composition of the (square) of the respective involutions isn't
        # perpendicular to the intersecting objects but contained in the plane.
        #
        return abs(self.half_trace.real()) < 1
    
    def is_not_stably_intersecting_or_touching_plane_at_inf(self):
        if sum(self.dimensions) < 3:
            return True

        return abs(self.half_trace.real()) > 1

    def dist(self):

        # Plane-plane: real trace
        # Plane-line: real trace (for square)
        # Plane-point: imag trace (sqrt(1 + square of imag part of half trace))
        # Line-line: complex trace 
        # Line-point: real trace (for square)
        # Point-point: real trace

        # Plane-plane and Plane-line:
        #    half trace < 1 means intersection
        #    half trace > 1 means neither touching nor intersecting
        #    otherwise: touching or identical

        t = self.half_trace

        if 0 in self.dimensions and 2 in self.dimensions:
            # Plane-point
            return arcsinh(abs(t.imag()))

        if self.dimensions[0] == 1 and self.dimensions[1] == 1:
            # Line-line

            # Use version of sqrt such that this evaluates even if t is close
            # to 1.
            return abs(log(abs(t + _branch_safe_sqrt(t * t - 1))))

        # Left:
        # Plane-plane
        # Plane-line
        # Line-point

        r = abs(t.real())
        RIF = r.parent()
        r = r.intersection(RIF(1, sage.all.Infinity))

        d = arccosh(r)

        if sum(self.dimensions) % 2 == 1:
            # We squared the composition matrix, so we need to
            # divide by two.
            return d / 2
        else:
            return d

    def fixed_projective_points_or_bounding_halfsphere(self):
        return fixed_projective_points_or_bounding_halfsphere_for_matrix(
            self.orientationPreservingComposition)

###############################################################################
# Various helpers

def _branch_safe_sqrt(z):
    if z != 0:
        return sqrt(z)

    d = sqrt(abs(z)).upper()
    RIF = z.real().parent()
    CIF = z.parent()

    return CIF(RIF(-d, d), RIF(-d, d))

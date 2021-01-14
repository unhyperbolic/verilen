from snappy.verify.upper_halfspace.finite_point import *
from .projectivePoint import *
from .triangleBase import *

__all__ = ['IdealTriangle']

class IdealTriangle(TriangleBase):
    """
    Represents an ideal triangle. Construct from three ProjectivePoint's.
    """

    def __init__(self, vertices):
        if not len(vertices) == 3:
            raise ValueError("Expected three vertices")
        for v in vertices:
            if not isinstance(v, ProjectivePoint):
                raise TypeError("Expected ProjectivePoint's instead of %r" % type(v))

        projective_endpoints = [
            TriangleBase.get_two(vertices, i) for i in range(3) ]
        super(IdealTriangle, self).__init__(
            vertices, projective_endpoints, [])

    def is_convex_hull_outside_extrusion(self, finitePoints):
        for finitePoint in finitePoints:
            if not isinstance(finitePoint, FinitePoint):
                raise TypeError("Expected FinitePoint")

        m = self._get_matrix_to_std_form()

        pts = [ pt.translate_PGL(m) for pt in finitePoints ]
        
        return (
            all([ pt.z.real() < 0 for pt in pts]) or
            all([ pt.z.real() > 1 for pt in pts]) or
            all([ (2 * pt.z.real() - 1) ** 2 + (2 * pt.t) ** 2 < 1 for pt in pts]))


            

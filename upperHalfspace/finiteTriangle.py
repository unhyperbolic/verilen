from snappy.verify.upper_halfspace.finite_point import *
from snappy.verify.upper_halfspace.ideal_point import _adjoint2
from .projectivePoint import *
from .reflections import *
from .triangleBase import *
from .primitivesDistanceEngine import *

import sage.all

__all__ = ['FiniteTriangle']

class FiniteTriangle(TriangleBase):
    """
    Represents a finite triangle. Construct from three FinitePoint's.
    """

    def __init__(self, vertices):
        if not len(vertices) == 3:
            raise ValueError("Expected three vertices")
        for v in vertices:
            if not isinstance(v, FinitePoint):
                raise TypeError("Expected FinitePoint's")

        vertex_reflections = [ 
            PointReflection.from_finite_point(v) for v in vertices ]

        projective_endpoints = [
            _get_projective_endpoints_of_line(
                TriangleBase.get_two(vertices, i),
                TriangleBase.get_two(vertex_reflections, i))
            for i in range(3) ]
        
        super(FiniteTriangle, self).__init__(
            vertices, projective_endpoints, vertex_reflections)

        self._matrices_to_and_vertices_in_std_plane = None

    def _get_matrices_to_and_vertices_in_std_plane(self):
        if self._matrices_to_and_vertices_in_std_plane is None:

            mstd = self._get_matrix_to_std_form()
            CIF = mstd.base_ring()

            matrices = [
                mstd,
                matrix(CIF, [[ 0, 1],
                             [-1, 1]]) * mstd,
                matrix(CIF, [[ 1,-1],
                             [ 1, 0]]) * mstd]

            self._matrices_to_and_vertices_in_std_plane = [
                (m, [ v.translate_PGL(m) for v in self.vertices ])
                for m in matrices ]
                
            for m, vertices in self._matrices_to_and_vertices_in_std_plane:
                for v in vertices:
                    if v.z.imag() != 0:
                        raise Exception(vertices)

        return self._matrices_to_and_vertices_in_std_plane

    def is_line_outside_extrusion(self, projectiveEndpoints):
        """
        Returns ``True`` if the line through the two given
        :class:`ProjectivePoint`\ s is verified to be outside the extrusion.

        Here, extrusion is the preimage of the finite triangle under
        orthogonal projection onto the plane containing the finite triangle.
        """

        if len(projectiveEndpoints) != 2:
            raise ValueError("Expected two ProjectivePoint's")

        for m, vertices in self._get_matrices_to_and_vertices_in_std_plane():
            r = _line_outside_finite_triangle(
                [ v.translate(m) for v in projectiveEndpoints ], vertices)
            if r > 0:
                return True
            if r < 0:
                return False

        return False


    """
    The old way of doing things...

        return self.are_vertices_outside_extrusion(idealVertices)

    def _are_vertices_outside_extrusion_one_side(self, vertices, i):
        j = (i+1) % 3
        for v in vertices:
            if not ProjectivePoint.cross_ratio(self.projective_endpoints[i][0],
                                               self.projective_endpoints[i][1],
                                               self.projective_endpoints[j][1],
                                               v).real() < 0:
                return False
        return True

    def are_vertices_outside_extrusion(self, vertices):
        for v in vertices:
            if not isinstance(v, ProjectivePoint):
                raise Exception("Expected projective point")

        for i in range(3):
            if self._are_vertices_outside_extrusion_one_side(vertices, i):
                return True

        return False
    """

###############################################################################
# Helpers

def _get_projective_endpoints_of_line(vertices, vertex_reflections):
    e = PrimitivesDistanceEngine([0,0], vertex_reflections)
    
    endpoints, h = e.fixed_projective_points_or_bounding_halfsphere()
    
    if not endpoints:
        raise Exception("Blah")

    CIF = vertices[0].z.parent()
    
    for i in range(1, 4):
        m = _adjoint2(
            ProjectivePoint.matrix_taking_0_1_inf_to_given_points(
                endpoints[0], ProjectivePoint(CIF(i)), endpoints[1]))
        if m.det().abs() > 0:
            try:
                vs = [ v.translate_PGL(m) for v in vertices ]
            except:
                continue
            
            if not all([not v.z != 0 for v in vs]):
                raise Exception("Images of vertices supposed to be "
                                "on vertical line.")

            if vs[0].t > vs[1].t:
                return endpoints[::-1]
            if vs[0].t < vs[1].t:
                return endpoints
    
    raise Exception("Could not find endpoints")

def _line_outside_finite_triangle(projectiveEndpoints, finitePoints):
    idealEndpoints = []
    for p in projectiveEndpoints:
        if p.inverted:
            if not p.z != 0:
                return 0
            idealEndpoints.append(1 / p.z)
        else:
            idealEndpoints.append(p.z)

    z = (idealEndpoints[0] + idealEndpoints[1]) / 2
    rSqr = _sqr(idealEndpoints[1] - idealEndpoints[0]) / 4

    signs = [ (_sqr(pt.z - z) + pt.t * pt.t) - rSqr for pt in finitePoints ]
    pos = [s > 0 for s in signs]
    neg = [s < 0 for s in signs]

    if all(pos) or all(neg):
        return +1

    if any(pos) and any(neg):
        return -1

    return 0

def _sqr(z):
    return z.real() * z.real() + z.imag() * z.imag()

################################################################################
#
# TESTING

class _FiniteTriangleTester(object):
    
    """
    A test rig for FiniteTriangle

    Run the test rig::
    
        sage: _FiniteTriangleTester().run_tests()

    """

    def run_tests(self):

        from sage.all import CIF, RIF
        
        pts = [
            FinitePoint(CIF(1.2,  5.2), RIF(2.4)),
            FinitePoint(CIF(3.2,  4.2), RIF(1.4)),
            FinitePoint(CIF(2.2, -3.2), RIF(0.4)) ]

        t = FiniteTriangle(pts)

        for i in range(3):
            m = t.reflections[1][i]
            if not abs(m[0,0] + m[1,1]) < 1e-6:
                raise Exception('%r' % m)

        for i, pt in enumerate(pts):
            if not pt.translate_PGL(t.reflections[0][i]).dist(pt) < 1e-6:
                raise Exception("Vertex reflection not fixing")

            for j in range(3):
                if i != j:
                    pt_t = pt.translate_PGL(t.reflections[1][j])
                    if not pt_t.dist(pt) < 1e-6:
                        raise Exception(
                            "Side reflection not fixing %r %r %r" % (
                                pt, pt_t, pt_t.dist(pt)))

            if not pt.translate_PGL(t.reflections[2][0]).dist(pt) < 1e-6:
                raise Exception("Plane reflection not fixing %r %r %r" % (
                        pt, pt.translate_PGL(t.reflections[2][0]),
                        pt.translate_PGL(t.reflections[2][0]).dist(pt)))

            

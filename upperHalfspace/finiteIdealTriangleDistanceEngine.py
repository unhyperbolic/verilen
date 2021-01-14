from .idealTriangle import *
from .finiteTriangle import *
from snappy.verify.upper_halfspace.finite_point import *
from .projectivePoint import *
from .fixedPoints import *
from .primitivesDistanceEngine import *

import itertools
import functools

__all__ = ['FiniteIdealTriangleDistanceEngine']

class FiniteIdealTriangleDistanceEngine(object):
    """
    Engine to compute a lower bound for the distance between a
    :class:`FiniteTriangle` and a :class:`IdealTriangle` that are assumed to
    not intersect.

    Example to compute lower bound for distance::

        sage: from sage.all import *
        sage: finiteTriangle = FiniteTriangle([
        ...        FinitePoint(CIF(0), RIF(2)),
        ...        FinitePoint(CIF(0.5), RIF(3)),
        ...        FinitePoint(CIF(-1), RIF(2.5))])
        sage: idealTriangle = IdealTriangle([
        ...        ProjectivePoint(CIF(0,1)),
        ...        ProjectivePoint(CIF(1,0)),
        ...        ProjectivePoint(CIF(-1,0))])
        sage: e = FiniteIdealTriangleDistanceEngine(finiteTriangle, idealTriangle)
        sage: e.compute_lowerbound()
        0.6931471805599453?

    """

    def __init__(self, finiteTriangle, idealTriangle):
        if not isinstance(finiteTriangle, FiniteTriangle):
            raise TypeError("Expected FiniteTriangle")
        if not isinstance(idealTriangle, IdealTriangle):
            raise TypeError("Expected IdealTriangle")

        self.triangles = [ finiteTriangle, idealTriangle ]
        self.finiteTriangle = finiteTriangle
        self.idealTriangle = idealTriangle

        self._distance_cells_cache = {}

    def compute_lowerbound(self):
        """
        Lower bound for the distance between the two triangles.
        """

        return self._compute_distance_cells(
            dimensions = (2, 2),
            cells = (0, 0),
            include_boundaries = (True, True))

    def _compute_distance_cells(self, dimensions, cells, include_boundaries):
        key = (dimensions, cells, include_boundaries)

        if not key in self._distance_cells_cache:
            self._distance_cells_cache[key] = self._compute_distance_cells_uncached(
                dimensions, cells, include_boundaries)

        return self._distance_cells_cache[key]

    def _compute_distance_cells_uncached(
        self, dimensions, cells, include_boundaries):

        d = self._compute_distance_primitives(dimensions, cells)
        if not d is None:
            return d

        boundaries = [
            triangle.get_boundary(dimension, cell) if include_boundary else []
            for triangle, dimension, cell, include_boundary
            in zip(self.triangles, dimensions, cells, include_boundaries) ]

        d = _min_allowing_none(
            [ self._compute_distance_cells(
                    _replace_in_tuple(dimensions, i, dimensions[i] - 1),
                    _replace_in_tuple(cells, i, boundary),
                    tuple([ i == j for j in range(2) ]))
              for i in range(2)
              for boundary in boundaries[i]])
        if not d is None:
            return d

        return _min_allowing_none(
            [ self._compute_distance_cells(
                    tuple([ dimension - 1 for dimension in dimensions ]),
                    boundary,
                    include_boundaries)
              for boundary in itertools.product(*boundaries) ])

    def _compute_distance_primitives(self, dimensions, cells):
        distanceEngine = PrimitivesDistanceEngine(
            dimensions,
            [ triangle.reflections[dimension][cell]
              for triangle, dimension, cell
              in zip(self.triangles, dimensions, cells) ])
        
        # We assume that the triangles are disjoint, so if the planes are
        # intersecting, then the shortest geodesic between them does not
        # hit the triangles in their interior. Thus report this as None.
        # Similar for plane and line.
        if distanceEngine.is_stably_intersecting():
            return None

        pts, b = distanceEngine.fixed_projective_points_or_bounding_halfsphere()
        if distanceEngine.is_not_stably_intersecting_or_touching_plane_at_inf():
            if pts:
                # We have two fixed points and they are the endpoints of the
                # shortest line segment between the two primitives
                if self._is_geodesic_outside_one_triangle(pts, dimensions, cells):
                    return None
                else:
                    return distanceEngine.dist()

        if sum(dimensions) == 3:
            # Line and plane can intersect perpendicular.
            # Detect this here since it is not handled by
            # distanceEngine.is_stably_intersecting
            lineIndex = dimensions.index(1)
            triangleIndex = dimensions.index(2)
            
            if FiniteIdealTriangleDistanceEngine._is_line_intersecting_plane_containing_triangle(
                    self.triangles[triangleIndex],
                    self.triangles[lineIndex].projective_endpoints[cells[lineIndex]]):
                return None

        if b:
            # The parabolic case
            for triangle, dimension, cell in zip(self.triangles, dimensions, cells):
                if b.are_outside(triangle.get_vertices_of_cell(dimension,cell)):
                    return None
            return None

        # The identity casse
        if dimensions[1] == 2:
            vs = self.finiteTriangle.get_vertices_of_cell(
                dimensions[0], cells[0])
            if self.idealTriangle.is_convex_hull_outside_extrusion(vs):
                return None
        elif dimensions[0] == 2:
            # Case plane - line
            vs = self.idealTriangle.get_vertices_of_cell(
                dimensions[1], cells[1])
            if self.finiteTriangle.is_line_outside_extrusion(vs):
                return None

        return distanceEngine.dist()

    @staticmethod
    def _is_line_intersecting_plane_containing_triangle(triangle, line_end_points):
        signs = [ 
            ProjectivePoint.cross_ratio(
                l, *triangle.three_projective_endpoints).imag()
            for l in line_end_points ]

        return signs[0] * signs[1] < 0

    def _is_geodesic_outside_one_triangle(self, pts, dimensions, cells):
        for triangle, dimension, cell in zip(self.triangles, dimensions, cells):

            boundary = triangle.get_boundary(dimension, cell)
            signs = triangle.get_sidedness_signs(pts, boundary)

            if any([s > 0 for s in signs]):
                return True

            if not all([s < 0 for s in signs]):
                if dimension == 1:
                    rpts = [ pt.translate(triangle.get_line_rotation(cell))
                             for pt in pts]
                    signs = triangle.get_sidedness_signs(rpts, boundary)
                    
                    if any([s > 0 for s in signs]):
                        return True

        return False

###############################################################################
# Various helpers

def _min(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return min(a, b)

def _min_allowing_none(values):
    if not values:
        return None
    return functools.reduce(_min, values)

def _replace_in_tuple(l, i, v):
    return tuple([ v if i == j else e for j, e in enumerate(l) ])

################################################################################
#
# TESTING

class _FiniteIdealTriangleDistanceEngineTester(object):
    """
    A test rig for FiniteIdealTriangleDistanceEngine.

    Run the test rig::
    
        sage: _FiniteIdealTriangleDistanceEngineTester().run_tests()

    """

    def example_plane_plane(self):
        from sage.all import exp

        d = self.RIF(1.125)

        i = IdealTriangle(
            [ ProjectivePoint(exp(self.CIF(d, 2.1 * i + 0.1)))
              for i in range(3)])

        e = FiniteIdealTriangleDistanceEngine(self.finiteTriangle1, i)
        
        if not abs(e.compute_lowerbound() - d) < 1e-20:
            raise Exception("Wrong distance")

    def example_plane_line(self):
        from sage.all import exp

        d = self.RIF(1.375)

        z = exp(self.CIF(d))

        for z1 in [ self.CIF(0), self.CIF(0.1, 0.1) ]:
            i = IdealTriangle(
                [ ProjectivePoint(z), ProjectivePoint(-z),
                  ProjectivePoint(z1, inverted = True) ])

            e = FiniteIdealTriangleDistanceEngine(self.finiteTriangle1, i)
        
            if not abs(e.compute_lowerbound() - d) < 1e-20:
                raise Exception("Wrong distance %r" % 
                                e.compute_lowerbound())

    def example_line_plane(self):
        from sage.all import exp, sin, cos

        d = self.RIF(1.09375)
        r = exp(d)

        alpha = self.RIF(0.4)
        z = self.CIF(r) * cos(alpha)

        for z1 in [ self.CIF(0), self.CIF(0.1, 0.1) ]:
            f = FiniteTriangle(
                [ FinitePoint(z =  z, t = r * sin(alpha)),
                  FinitePoint(z = -z, t = r * sin(alpha)),
                  FinitePoint(z1, t = self.RIF(10.0)) ])

            e = FiniteIdealTriangleDistanceEngine(f, self.idealTriangle1)
        
            if not abs(e.compute_lowerbound() - d) < 1e-20:
                raise Exception("Wrong distance %r" % 
                                e.compute_lowerbound())

    def finitePointFromPolar(self, d, alpha, beta):
        from sage.all import exp, sin, cos

        return FinitePoint(
            z = sin(alpha) * exp(self.CIF(d, beta)),
            t = cos(alpha) * exp(d))

    def example_line_line(self):
        from sage.all import exp, sin, cos

        d = self.RIF(1.09375)
        r = exp(d)

        alpha = self.RIF(0.4)
        
        for l in [ 0, 0.8, 1.2, self.RIF.pi() / 2]:
            z = self.CIF(r) * exp(self.CIF(0,l)) * cos(alpha)
            
            for z1 in [ self.CIF(0), self.CIF(0.1, 0.1) ]:

                f = FiniteTriangle(
                    [ FinitePoint(z =  z, t = r * sin(alpha)),
                      FinitePoint(z = -z, t = r * sin(alpha)),
                      FinitePoint(z1, t = self.RIF(10.0)) ])
                
                e = FiniteIdealTriangleDistanceEngine(f, self.idealTriangle2)
                
                if not abs(e.compute_lowerbound() - d) < 1e-20:
                    raise Exception("Wrong distance %r" % 
                                    e.compute_lowerbound())

    def example_point_plane(self):
        d = self.RIF(1.75)

        f = FiniteTriangle(
            [ self.finitePointFromPolar(d, self.RIF(0), self.RIF(0)),
              self.finitePointFromPolar(2 * d, self.RIF(0.1), self.RIF(1)),
              self.finitePointFromPolar(2 * d, self.RIF(0.1), self.RIF(2))])
        
        e = FiniteIdealTriangleDistanceEngine(f, self.idealTriangle1)

        if not abs(e.compute_lowerbound() - d) < 1e-20:
            raise Exception("Wrong distance %r" % 
                            e.compute_lowerbound())

    def example_point_line(self):
        d = self.RIF(1.625)

        f = FiniteTriangle(
            [ self.finitePointFromPolar(2 * d, self.RIF(0), self.RIF(0)),
              self.finitePointFromPolar(d,  self.RIF(0.1), self.RIF(1)),
              self.finitePointFromPolar(d, -self.RIF(0.1), self.RIF(1))])
        
        e = FiniteIdealTriangleDistanceEngine(f, self.idealTriangle1)

        if not abs(e.compute_lowerbound() - d) < 1e-20:
            raise Exception("Wrong distance %r" % 
                            e.compute_lowerbound())

    def example_point_plane_parallel(self):
        a = self.RIF(1.1875)

        ref_pt = self.finitePointFromPolar(a, self.RIF(0), self.RIF(0))
        pt0 = self.finitePointFromPolar(a, self.RIF(0.3), self.RIF.pi()/2)
        pt1 = FinitePoint(self.CIF( 0.4, pt0.z.imag()), self.RIF(0.5))
        pt2 = FinitePoint(self.CIF(-0.4, pt0.z.imag()), self.RIF(0.5))
        f = FiniteTriangle([pt0, pt1, pt2])

        d = ref_pt.dist(pt0)
        
        e = FiniteIdealTriangleDistanceEngine(f, self.idealTriangle3)

        if not abs(e.compute_lowerbound() - d) < 1e-20:
            raise Exception("Wrong distance %r" % 
                            e.compute_lowerbound())
        
    def example_t(self):
        from sage.all import exp, sin, cos

        d = self.RIF(1.09375)

        f = FiniteTriangle(
            [ self.finitePointFromPolar(d, 0.4, 0.2),
              self.finitePointFromPolar(d, 0.6, 0.2),
              self.finitePointFromPolar(d, 0.5, 0.3) ])
        
        e = FiniteIdealTriangleDistanceEngine(f, self.idealTriangle2)
        
        if not abs(e.compute_lowerbound() - d) < 1e-20:
            raise Exception("Wrong distance %r" % 
                            e.compute_lowerbound())


    def setup_fields_and_triangles(self):
        from sage.all import RealIntervalField, ComplexIntervalField, sin, cos, exp

        self.RIF = RealIntervalField(120)
        self.CIF = ComplexIntervalField(120)

        alpha = self.RIF(0.4)

        self.finiteTriangle1 = FiniteTriangle(
            [ FinitePoint(z = exp(self.CIF(0, 2.1 * i + 0.1)) * cos(alpha),
                          t = sin(alpha))
              for i in range(3) ])

        self.idealTriangle1 = IdealTriangle(
            [ ProjectivePoint(z = exp(self.CIF(0, 2.1 * i + 0.1)))
              for i in range(3) ])

        self.idealTriangle2 = IdealTriangle(
            [ ProjectivePoint(self.CIF(z)) for z in [-1, 0, 1]])

        self.idealTriangle3 = IdealTriangle(
            [ ProjectivePoint(self.CIF( 1)),
              ProjectivePoint(self.CIF(-1)),
              ProjectivePoint(self.CIF(0), True)])

    def my_example(self):
        from sage.all import CIF, RIF
        finiteTriangle = FiniteTriangle([
                FinitePoint(CIF(0), RIF(2)),
                FinitePoint(CIF(0.5), RIF(3)),
                FinitePoint(CIF(-1), RIF(2.5))])
        idealTriangle = IdealTriangle([
                ProjectivePoint(CIF(0,1)),
                ProjectivePoint(CIF(1,0)),
                ProjectivePoint(CIF(-1,0))])

        #print finiteTriangle.reflections[2][0]
        #print idealTriangle.reflections[1][0]

        e = FiniteIdealTriangleDistanceEngine(finiteTriangle, idealTriangle)
        d = e.compute_lowerbound()

        if d < 0.0001:
            raise Exception()
        
    def run_tests(self):
        self.setup_fields_and_triangles()

        self.example_plane_plane()
        self.example_plane_line()
        self.example_line_plane()
        self.example_line_line()
        self.example_point_plane()
        self.example_point_line()
        self.example_point_plane_parallel()

        self.my_example()


from .reflections import *

from snappy.verify.upper_halfspace.ideal_point import _adjoint2

class TriangleBase(object):
    """
    Base class for :class:`FiniteTriangle` and :class:`IdealTriangle`. Some of
    the data stored:

       * ``vertices`` are the finite or ideal vertices of the triangle
       * ``projective_endpoints`` contains for each edge the two
         endpoints (of type :class:`ProjectivePoint`) of the line containing
         the respective edge.
       * ``three_projective_endpoints`` contains three of the (not necessarily)
         distinct six endpoints. They span the plane containing the triangle.
       * ``reflections`` stores for each dimension and for each cell of that
         dimension the involution fixing that cell
    """

    def __init__(self, vertices, projective_endpoints, vertex_reflections):
        self.vertices = vertices
        self.projective_endpoints = projective_endpoints

        self.three_projective_endpoints = [
            pts[0] for pts in projective_endpoints ]

        self.reflections = [
            vertex_reflections, 
            [ LineReflection.from_two_projective_points(*pts)
              for pts in projective_endpoints ],
            [ PlaneReflection.from_three_projective_points(
                    *self.three_projective_endpoints) ] ]

        # Data we cache
        self._line_rotations = [ None, None, None ]
        self._matrix_to_std_form = None

    def _get_matrix_to_std_form(self):
        """
        Compute the matrix bringing the triangle into the plane with imaginary
        part 0.

        For an ideal triangle, the matrix will bring the vertices to 0, 1, and
        infinity.
        """

        if self._matrix_to_std_form is None:
            self._matrix_to_std_form = _adjoint2(
                ProjectivePoint.matrix_taking_0_1_inf_to_given_points(
                    *self.three_projective_endpoints))
        return self._matrix_to_std_form

    def _orient_line(self, pts):
        """
        Since the vertices of triangle have a cyclic ordering, the plane
        containing the triangle naturally has a "positive" and a "negative"
        side from which the plane appears to be oriented positively,
        respectively, negatively. 
        
        Given two :class:`ProjectivePoint`\ s, flip them if necessary so that
        the second one is always to the same side.

        Note that it is not checked that the first point is to the other side
        of the plane.
        """

        im = ProjectivePoint.cross_ratio(
            pts[1], *self.three_projective_endpoints).imag()

        if im > 0:
            return pts
        if im < 0:
            return pts[::-1]

        return None

    def get_sidedness_signs(self, pts, sides):
        opts = self._orient_line(pts)
        if not opts:
            RIF = pts[0].z.real().parent()
            return [ RIF(-1, 1) for i in sides ]

        return [
            ProjectivePoint.cross_ratio(*(opts + self.projective_endpoints[i])).imag()
            for i in sides ]

    @staticmethod
    def get_two(l, i):
        return [ l[(i + 1) % 3], l[(i + 2) % 3] ]

    def get_boundary(self, dimension, cell):
        if dimension == 2:
            return range(3)
        if dimension == 1 and self.reflections[0]:
            return TriangleBase.get_two(list(range(3)), cell)
        return []

    def get_vertices_of_cell(self, dimension, cell):
        if dimension == 2:
            return self.vertices
        if dimension == 1:
            return TriangleBase.get_two(self.vertices, cell)
        if dimension == 0:
            return [ self.vertices[cell] ]

    def get_line_rotation(self, i):
        if self._line_rotations[i] is None:
            self._line_rotations[i] = LineRotation.from_two_projective_points(
                *self.projective_endpoints[i])
        return self._line_rotations[i]


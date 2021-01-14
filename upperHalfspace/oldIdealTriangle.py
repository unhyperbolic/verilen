from .projectivePoint import *
from .idealPoint import _adjoint2

__all__ = ['IdealTriangle']

class IdealTriangle(object):
    def __init__(self,
                 baseIdealTriangle = None,
                 matrix = None,
                 side_reflections = 3 * [ None ],
                 plane_reflection = None,
                 matrix_to_std_form = None):
        self._baseIdealTriangle = baseIdealTriangle
        self._matrix = matrix
        self._side_reflections = side_reflections
        self._plane_reflection = plane_reflection
        self._matrix_to_std_form = matrix_to_std_form

    @staticmethod
    def from_vertices(vertices):
        return IdealTriangle(
            side_reflections = [
                ProjectivePoint.line_reflection_fixing_points(
                    vertices[ (i+1) % 3], vertices[ (i+2) % 3])
                for i in range(3) ],
            plane_reflection = (
                ProjectivePoint.plane_reflection_fixing_points(
                    *vertices)),
            matrix_to_std_form = _adjoint2(
                ProjectivePoint.matrix_taking_0_1_inf_to_given_points(
                    *vertices)))

    def translate(self, matrix):
        return IdealTriangle(
            baseIdealTriangle = self,
            matrix = matrix)

    def get_side_reflection(self, i):
        if self._side_reflections[i] is None:
            self._side_reflections[i] = (
                matrix *
                baseIdealTriangle.get_side_reflection(i) *
                _adjoint2(matrix))
        return self._side_reflections[i]

    def get_plane_reflection(self):
        if self._plane_reflection is None:
            self._plane_reflection = (
                matrix *
                baseIdealTriangle.get_plane_reflection() *
                _adjoint2(matrix))
        return self._plane_reflection

    def get_matrix_to_std_form(self):
        if self._matrix_to_std_form is None:
            self._matrix_to_std_form = (
                baseIdealTriangle.get_matrix_to_std_form() *
                _adjoint2(matrix))
        return self._matrix_to_std_form

    def signed_distance_measure(self, finitePoint):
        pt = finitePoint.translate_PGL(self.get_matrix_to_std_form())
        
        return interval_aware_max(
            [ pt.z.real() - 1,
              -pt.z.real(),
              1 - sqrt((2 * pt.z.real() - 1) ** 2 + (2 * pt.t) ** 2) ])

    def is_point_outside(self, finitePoint):
        return self.signed_distance_measure(finitePoint) > 0

    def finite_triangle_outside_extrusion(self, finitePoints):
        pts = [ pt.translate_PGL(self.get_matrix_to_std_form())
                for pt in finitePoints ]
        
        return (all([ pt.z.real() < 0 for pt in pts]) or
                all([ pt.z.real() > 1 for pt in pts]) or
                all([ (2 * pt.z.real() - 1) ** 2 + (2 * pt.t) ** 2 < 1 for pt in pts]))

from .projectivePoint import *
from snappy.verify.upper_halfspace.ideal_point import *
from snappy.verify.upper_halfspace.finite_point import *
from snappy.verify.upper_halfspace.extended_matrix import *
from .reflections import *
from .boundingHalfsphere import *

import sage.all
from sage.all import matrix, sqrt

__all__ = [
    'fixed_projective_points_or_bounding_halfsphere_for_matrix'
    ]

def fixed_projective_points_or_bounding_halfsphere_for_matrix(m):
    """
    Fixed points

        sage: from sage.all import *
        sage: m = matrix([[CIF(2.3, 3.1), CIF(1.2,-0.1)],
        ...               [CIF(5.2, 1.2), CIF(6.2, 1.2)]])
        sage: fixed_projective_points_or_bounding_halfsphere_for_matrix(m) # doctest: +NUMERIC12
        ([ProjectivePoint(-0.858660982188960? + 0.467887984415522?*I), ProjectivePoint(0.226638510278847? + 0.043347970640658?*I)], BoundingHalfsphere(ProjectivePoint(-0.316011235955057? + 0.2556179775280899?*I), 0.582689714094655))
        
    """

    if isinstance(m, ExtendedMatrix):
        if m.isOrientationReversing:
            raise ValueError(
                "This function does not support computing fixed "
                "points of orientation-reversing extended matrices.")
        m = m.matrix

    if m[1,0] != 0:
        return _fixed_ideal_points_for_matrix(m, inverted = False)

    if m[0,1] != 0:
        m = matrix([[m[1,1],m[1,0]],[m[0,1],m[0,0]]])
        return _fixed_ideal_points_for_matrix(m, inverted = True)

    h    = matrix(m.base_ring(), [[ 1, 0], [ 1, 1]])
    hinv = matrix(m.base_ring(), [[ 1, 0], [-1, 1]])

    m = hinv * m * h

    if m[1, 0] != 0:
        pts, b = _fixed_ideal_points_for_matrix(m, inverted = False)
        if pts:
            return [ pt.translate(h) for pt in pts], None

    return None, None

################################################################################
#
# Various helper functions
    
def _fixed_ideal_points_for_matrix(m, inverted):
    # (a * z + b) = z * (c * z + d)
    # z^2 + (d - a) / c * z - b / c = 0

    # p' = (a - d) / 2 * c
    # q = - b / c

    # p' +/- sqrt( p'^2 - q ) = -p' +/- sqrt(p'^2 + b / c)

    cinv = 1 / m[1, 0]

    p = (m[0,0] - m[1,1]) * cinv / 2
    d = p ** 2 + m[0,1] * cinv

    r = d.abs().sqrt().upper()
    b = BoundingHalfsphere(ProjectivePoint(p, inverted), r)

    if d != 0:
        s = d.sqrt()
        result = [p - s, p + s]
        if result[0] != result[1]:
            return [ ProjectivePoint(z, inverted) for z in result ], b

    return None, b

################################################################################
#
# TESTING

class _FixedPointsTester(object):

    """
    A test rig for fixedPoints

    Run the test rig::

        sage: _FixedPointsTester().run_tests()

    """

    def matrices(self):
        from sage.all import RealIntervalField, ComplexIntervalField, matrix
        RIF = RealIntervalField(120)
        CIF = ComplexIntervalField(120)

        return [
            matrix(
                [[CIF(RIF(1.3),RIF(-0.4)), CIF(RIF(5.6),RIF(2.3))],
                 [CIF(RIF(-0.3), RIF(0.1)), CIF(1)]]),
            matrix(
                [[CIF(RIF(0.3),RIF(-1.4)), CIF(RIF(3.6),RIF(6.3))],
                 [CIF(RIF(-0.3), RIF(1.1)), CIF(1)]]),
            matrix(
                [[CIF(2.3, 1.2), CIF(0)],
                 [CIF(0), CIF(1)]]),
            matrix(
                [[CIF(2.3, 1.2), CIF(5)],
                 [CIF(0), CIF(1)]]),
            matrix(
                [[CIF(2.3, 1.2), CIF(0)],
                 [CIF(5.6), CIF(1)]]),
            matrix(
                [[CIF(1), CIF(10.0)],
                 [CIF(0), CIF(1)]]),
            matrix(
                [[CIF(1), CIF(0)],
                 [CIF(10.0), CIF(1)]]),
            matrix(
                [[CIF(0), CIF(-1)],
                 [CIF(1), CIF(2)]]),
            matrix(
                [[CIF(2), CIF(-1)],
                 [CIF(1), CIF(0)]])
            ]

    def run_tests(self):
        self.test_fixed_points_for_matrix()
        
    def test_fixed_points_for_matrix(self):
        for m in self.matrices():
            pts, h = fixed_projective_points_or_bounding_halfsphere_for_matrix(m)
            if not m[0,0] + m[1,1] != 2:
                pts = [ h.projectivePoint]
            for pt in pts:
                if not pt.translate(m).is_close_to(pt):
                    raise Exception("Point not fixed %r %r" % (
                            pt.translate(m), pt))
                    


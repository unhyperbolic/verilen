from snappy.verify.upper_halfspace.finite_point import *
from .projectivePoint import *

import sage.all
from sage.all import matrix

__all__ = ['BoundingHalfsphere']

class BoundingHalfsphere:
    """
    A halfspace in hyperbolic space represented by a :class:`ProjectivePoint`
    ``p`` on the boundary of hyperbolic space and a Euclidean radius ``r``
    (expected to be an element in SageMath's ``RealIntervalField``).

    Take the halfspace cut out by the hyperbolic halfplane corresponding to
    a halfsphere with the given Euclidean radius about ``p.z`` in the upper half
    plane model. If the ``inverted`` flag is set to ``False``, this is the 
    halfspace we represent here, otherwise, let 1/z act on it by Poincare
    extension.
    """

    def __init__(self, projectivePoint, r):
        self.projectivePoint = projectivePoint
        self.r = r
        self._rSqr = r * r

    def are_outside(self, pts):
        """
        Given a list of :class:`FinitePoint`\ s or :class:`ProjectivePoint`\ s,
        returns ``True`` if each is outside the halfspace::

            sage: from sage.all import *
            sage: b = BoundingHalfsphere(ProjectivePoint(CIF(2,1)), RIF(1))
            sage: b.are_outside([FinitePoint(CIF(3,2),RIF(0.5)), FinitePoint(CIF(4,5), RIF(2))])
            True
            sage: b.are_outside([ProjectivePoint(CIF(1.9, 0.9))])
            False

        """
        for pt in pts:
            if not self.is_outside(pt):
                return False
        return True

    def is_outside(self, pt):
        """
        Given a single :class:`FinitePoint` or :class:`ProjectivePoint`, returns
        ``True`` if it is outside the halfspace::


            sage: from sage.all import *
            sage: b = BoundingHalfsphere(ProjectivePoint(CIF(2,0), inverted = True), RIF(1))
            sage: b.is_outside(FinitePoint(CIF(3,2),RIF(0.5)))
            True
            sage: b.is_outside(ProjectivePoint(CIF(1.9, 0.9), inverted = True))
            False
            sage: b.is_outside(ProjectivePoint(CIF(0.45, 0.1)))
            False
        """

        if isinstance(pt, FinitePoint):
            if self.projectivePoint.inverted:
                pt = pt.translate_PSL(self._get_inversion_matrix())

            dSqr = _sqr(pt.z - self.projectivePoint.z) + pt.t * pt.t

        elif isinstance(pt, ProjectivePoint):
            if self.projectivePoint.inverted ^ pt.inverted:
                if not pt.z != 0:
                    return False
                z = 1 / pt.z
            else:
                z = pt.z

            dSqr = _sqr(z - self.projectivePoint.z)
        else:
            raise TypeError(
                "Expected FinitePoint or ProjectivePoint, not %r" % pt)

        return dSqr > self._rSqr

    def __repr__(self):
        return "BoundingHalfsphere(%r, %r)" % (self.projectivePoint, self.r)

    def _get_inversion_matrix(self):
        CIF = self.projectivePoint.z.parent()
        return matrix(CIF, [[ 0, sage.all.I],
                            [ sage.all.I, 0]])

###############################################################################
# Various helpers

def _sqr(z):
    return z.real() * z.real() + z.imag() * z.imag()

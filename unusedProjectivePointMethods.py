from upperHalfspace import *
from upperHalfspace import _numerator_of_diff

def fourth_point(z, z0, z1, z2):
    """
    Given a cross ratio z and three ProjectivePoint's z0, z1, z2, compute
    the fourth ProjectivePoint z3 such that (z0, z1, z2, z3) have cross ratio z.
    z is supposed to be an element in ``ComplexIntervalField``.

    >>> from sage.rings.real_mpfi import RealIntervalField
    >>> from sage.rings.complex_interval_field import ComplexIntervalField
    >>> from sage.matrix.constructor import MatrixFactory
    >>> RIF = RealIntervalField()
    >>> CIF = ComplexIntervalField()
    >>> matrix = MatrixFactory()
    
    >>> a = CIF(RIF(1.3,1.30000000001),RIF(-0.4,-0.3999999999999))
    >>> b = CIF(RIF(5.6,5.60000000001),RIF(2.3,2.30000000001))
    >>> c = CIF(RIF(-0.3,-0.29999999999), RIF(0.1,0.10000000001))
    >>> d = (1 + b * c) / a

    >>> m = matrix([[a,b], [c,d]])

    >>> z = CIF(1.2, 1.3)
    >>> z0 = ProjectivePoint(CIF(0))
    >>> z1 = ProjectivePoint(CIF(0), inverted = True)
    >>> z2 = ProjectivePoint(CIF(1))

    >>> for i in range(6):
    ...     print ProjectivePoint.cross_ratio(z0, z1, z2, fourth_point(z, z0, z1, z2))
    ...     z0 = z0.translate(m)
    ...     z1 = z1.translate(m)
    ...     z2 = z2.translate(m)
    1.2000000000000000? + 1.3000000000000001?*I
    1.2000000? + 1.3000000?*I
    1.200000? + 1.300000?*I
    1.20000? + 1.30000?*I
    1.2000? + 1.3000?*I
    1.200? + 1.300?*I
    
    """
    
    # Solving the above equation for z3 yields
    #
    #       z * (z  - z ) * z  - (z  - z ) * z
    #             2    1     0     2    0     1
    # z  = -------------------------------------
    #  3    z * (z  - z )      - (z  - z )
    #             2    1           2    0
    #
    # This again simplifies after substituiting to:
    
    return ProjectivePoint.fromNumeratorAndDenominator(
          z * _numerator_of_diff(z2, z1) * z0.get_numerator()
            - _numerator_of_diff(z2, z0) * z1.get_numerator(),
          z * _numerator_of_diff(z2, z1) * z0.get_denominator()
            - _numerator_of_diff(z2, z0) * z1.get_denominator())

def fill_missing_vertex_for_tet(z, verts):
    """
    Given a list ``verts`` of four elements, three of them
    :class:`ProjectivePoint` and one ``None``, replace ``None`` in the list
    by an :class:`ProjectivePoint` such that their cross ratio becomes the given
    z (an element in ``ComplexIntervalField``).
    The function has no return value as ``verts`` is modified in place.

    >>> from sage.all import *
    >>> z = CIF(0.31,0.5)
    >>> verts = [ProjectivePoint(CIF(0)), ProjectivePoint(CIF(0),inverted=True), ProjectivePoint(CIF(1)), None]
    >>> fill_missing_vertex_for_tet(z, verts)
    >>> verts
    [ProjectivePoint(0), ProjectivePoint(0, inverted = True), ProjectivePoint(1), ProjectivePoint(0.31000000000000000? + 0.50000000000000000?*I, inverted = True)]

    >>> verts = [ProjectivePoint(CIF(1.2,3.4)), ProjectivePoint(CIF(3.4,1.2)), None, ProjectivePoint(CIF(6.2,1.2))]
    >>> fill_missing_vertex_for_tet(z, verts)
    >>> ProjectivePoint.cross_ratio(*verts)
    0.31000000000000? + 0.50000000000000?*I

    >>> verts = [ProjectivePoint(CIF(1.2,3.4)), None, ProjectivePoint(CIF(3.4,1.2)), ProjectivePoint(CIF(6.2,1.2))]
    >>> fill_missing_vertex_for_tet(z, verts)
    >>> ProjectivePoint.cross_ratio(*verts)
    0.31000000000000? + 0.50000000000000?*I
    
    >>> verts = [None, ProjectivePoint(CIF(3.2,1.2)), ProjectivePoint(CIF(3.4,1.2)), ProjectivePoint(CIF(6.2,1.2))]
    >>> fill_missing_vertex_for_tet(z, verts)
    >>> ProjectivePoint.cross_ratio(*verts)
    0.3100000000000? + 0.50000000000000?*I
    """

    # By applying a proper permutation to the vertices of a simplex,
    # we can assume that the missing vertex is z3 and apply the above
    # method to compute it.
    
    # Such a permutation in S_4 should not change the cross ratio, i.e.,
    # it should act trivially on pairs of opposite edges.
    # The kernel ker(S_4->S_3) of the action on pairs of opposite edges
    # is the Klein four group Z/2 x Z/2.

    if verts[0] is None:
        # (03)(12)
        verts[0] = fourth_point(z, verts[3], verts[2], verts[1])
    if verts[1] is None:
        # (02)(13)
        verts[1] = fourth_point(z, verts[2], verts[3], verts[0])
    if verts[2] is None:
        # (01)(23)
        verts[2] = fourth_point(z, verts[1], verts[0], verts[3])
    if verts[3] is None:
        # id
        verts[3] = fourth_point(z, verts[0], verts[1], verts[2])

def _to_special_form(m):
    for i in [0, 1]:
        max_index, max_val = max(enumerate(m.column(i)[i:]),
                                 key = lambda x:x[1].abs().lower())

        if max_val.contains_zero():
            raise ZeroDivisionError

        if max_index != 0:
            m[max_index + i], m[i] = m[i], m[max_index + i]

        m[i] /= m[i][i]

        for j in range(i+1, 3):
            m[j] -= m[j][i] * m[i]

def matrix_taking_triple_to_triple(srcProjectivePoints, dstProjectivePoints):
    
    CIF = srcProjectivePoints[0].z.parent()

    eqns = matrix(
        CIF,
        [ [   srcPt.get_numerator()   * dstPt.get_denominator(),
              srcPt.get_denominator() * dstPt.get_denominator(),
            - srcPt.get_numerator()   * dstPt.get_numerator(),
            - srcPt.get_denominator() * dstPt.get_numerator() ]
          for srcPt, dstPt in zip(srcProjectivePoints, dstProjectivePoints) ])
    
    _to_special_form(eqns)

    c =                                    -     eqns[2][3]
    d =                         eqns[2][2]
    b = - (                 c * eqns[1][2] + d * eqns[1][3])
    a = - (b * eqns[0][1] + c * eqns[0][2] + d * eqns[0][3])
    
    return matrix(CIF, [[a,b],[c,d]])

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()


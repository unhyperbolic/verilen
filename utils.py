
from sage.all import sqrt

def tangent_vector_finite_points(finPoint1, finPoint2):
    
    diffZ = finPoint2.z - finPoint1.z
    absSqrDiffZ = diffZ.real() ** 2 + diffZ.imag() ** 2

    n = (finPoint1.t ** 2 - finPoint2.t ** 2 - absSqrDiffZ) / 2
    
    vLen = sqrt(absSqrDiffZ + (n / finPoint1.t) ** 2)

    return finPoint1.t * diffZ / vLen, -n / vLen

    
    

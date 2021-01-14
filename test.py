import upperHalfspace.projectivePoint
import upperHalfspace.fixedPoints
import upperHalfspace.finiteTriangle
import upperHalfspace.reflections
import upperHalfspace.boundingHalfsphere
import upperHalfspace.finiteIdealTriangleDistanceEngine
import upperHalfspace.primitivesDistanceEngine

import simpleTilingEngine

from snappy.sage_helper import doctest_modules

import sys, getopt

modules = [
    upperHalfspace.projectivePoint,
    upperHalfspace.fixedPoints,
    upperHalfspace.finiteTriangle,
    upperHalfspace.reflections,
    upperHalfspace.boundingHalfsphere,
    upperHalfspace.finiteIdealTriangleDistanceEngine,
    upperHalfspace.primitivesDistanceEngine,
    simpleTilingEngine,
    ]


def run_doctests(verbose=False, print_info=True):
    return doctest_modules(modules,
                           verbose = verbose, print_info = print_info)

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    run_doctests(verbose)

import sage.all
from sage.all import sqrt
from upperHalfspace import *
from snappy.verify.upper_halfspace.ideal_point import *
from upperHalfspace.projectivePoint import *
from snappy.snap.t3mlite import simplex

def matrices_are_close(m1, m2, epsilon = 1e-6):
    m1 = m1 / sqrt(m1.det())
    m2 = m2 / sqrt(m2.det())

    return any(
        [ all([abs(m1[i,j] - s * m2[i,j]) < 1e-6
               for i in range(2) for j in range(2)])
          for s in [+1, -1]])

def simple_tiling_engine_consistency_check(simpleTilingEngine):

    tiling_engine_base_consistency_check(simpleTilingEngine)

    tet = simpleTilingEngine.mcomplex.Tetrahedra[0]
    CIF = tet.ShapeParameters[simplex.E01].parent()
    RIF = tet.ShapeParameters[simplex.E01].real().parent()
        
    tiles = simpleTilingEngine.intervalTree.find(
        RIF(-sage.all.Infinity, sage.all.Infinity))

    for tile in tiles:
        for tet in simpleTilingEngine.mcomplex.Tetrahedra:
            idealPoints = [
                tile.vertexToIdealPoint[tet.Class[v]]
                for v in simplex.ZeroSubsimplices ]

            z = cross_ratio(*idealPoints)

            if not (abs(z - tet.ShapeParameters[simplex.E01]) < 1e-6):
                print()
                print(idealPoints)
                print(z)
                print(tet.ShapeParameters[simplex.E01])

        for vertex in simpleTilingEngine.mcomplex.Vertices:
            tilingVertex = tile.vertexToIdealPoint[vertex]

            idealPoint = ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(
                CIF, vertex.IdealPoint)

            idealPointImg = idealPoint.translate(tile.matrix)
                
            tilingIdealPoint = ProjectivePoint.fromComplexIntervalFieldAndIdealPoint(
                CIF, tilingVertex)

            if not idealPointImg.is_close_to(tilingIdealPoint):
                print(idealPointImg, tilingIdealPoint)

def tiling_engine_base_consistency_check(tilingEngineBase):

    tet = tilingEngineBase.mcomplex.Tetrahedra[0]
    CIF = tet.ShapeParameters[simplex.E01].parent()
    RIF = tet.ShapeParameters[simplex.E01].real().parent()
        
    tiles = tilingEngineBase.intervalTree.find(
        RIF(-sage.all.Infinity, sage.all.Infinity))

    for tile in tiles:

        for g, m in tilingEngineBase.mcomplex.GeneratorMatrices.items():
            other_tile = tilingEngineBase.find_tile(tile.matrix * m)
                
            neighboringTile = tile.generatorToNeighboringTile.get(g)

            if not neighboringTile == other_tile:
                print("BAD")
                
            if neighboringTile:
                if not matrices_are_close(
                    neighboringTile.matrix, tile.matrix * m):
                    print("BADBAD")
                                           

from snappy.snap.mcomplex_base import *

from snappy.snap.fundamental_polyhedron import *
from snappy.snap.kernel_structures import *
from snappy.snap.t3mlite import simplex
from snappy.snap import t3mlite as t3m

from snappy.verify.upper_halfspace.ideal_point import *
from upperHalfspace.finiteTriangle import *

from sage.all import prod, Rational

__all__ = ['SpineEngine']

class SpineEngine(McomplexEngine):
    """

    >>> from snappy import *
    >>> M = Manifold("o9_24124")
    >>> shapes = M.verify_hyperbolicity(bits_prec = 512)[1]
    >>> t = SpineEngine.from_manifold_and_shapes(M, shapes)
    >>> t._check_consistency(epsilon = 1e-50)

    """

    class EdgeInfo(object):
        def __init__(self, perm, offset):
            self.perm = perm
            self.offset = offset

    @staticmethod
    def from_manifold_and_shapes(
            manifold, shapes, normalize_matrices = False, match_kernel = True):
        m = t3m.Mcomplex(manifold)
        
        t = TransferKernelStructuresEngine(m, manifold)
        t.add_shapes(shapes)
        t.choose_and_transfer_generators(
            compute_corners = True, centroid_at_origin = False)
        
        s = SpineEngine(m)
        s.compute_edge_info_before_ungluing()
        
        f = FundamentalPolyhedronEngine(m)
        f.unglue()

        if match_kernel:
            init_verts = f.init_vertices_kernel()
        else:
            init_verts = f.init_vertices()

        f.visit_tetrahedra_to_compute_vertices(
            m.ChooseGenInitialTet, init_verts)
        f.compute_matrices(normalize_matrices = normalize_matrices)

        s.add_spine_after_ungluing()
        s.add_finite_spine_triangles()

        return s

    def __init__(self, mcomplex):
        super(SpineEngine, self).__init__(mcomplex)

        self._edgeInfos = [ {} for tet in mcomplex.Tetrahedra ]

    def compute_edge_info_before_ungluing(self):
        for edge in self.mcomplex.Edges:
            self._compute_edge_info_for_edge(edge)

    def _RIF(self):
        tet = self.mcomplex.Tetrahedra[0]
        shape = tet.ShapeParameters[simplex.E01]
        return shape.real().parent() 

    def _compute_edge_info_for_edge(self, edge):
        RIF = self._RIF()
        distance_of_vertex = RIF(1)
        distances_of_vertices = []

        for tet, perm in edge.embeddings():
            distances_of_vertices.append(distance_of_vertex)

            z = tet.ShapeParameters[perm.image(simplex.E01)]
            distance_of_vertex *= abs(z)

        offsets = _invert_and_normalize(distances_of_vertices)
            
        for offset, (tet, perm) in zip(offsets, edge.embeddings()):
            self._edgeInfos[tet.Index][perm.image(simplex.E01)] = (
                SpineEngine.EdgeInfo(perm, offset))
    
    def add_spine_after_ungluing(self):
        self._add_tet_centers()
        self._add_face_centers()
        self._add_edge_centers()

    def add_finite_spine_triangles(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.spineTriangles = {
                perm : FiniteTriangle([
                        tet.InCenter,
                        tet.Class[perm.image(simplex.F3 )].InCenter,
                        tet.Class[perm.image(simplex.E01)].MidPoint])
                for perm in t3m.Perm4.A4() }
        
    def _add_tet_centers(self):
        for tet in self.mcomplex.Tetrahedra:
            tet.InRadius, tet.InCenter = compute_inradius_and_incenter(
                [ tet.Class[v].IdealPoint for v in simplex.ZeroSubsimplices ])
            
    def _add_face_centers(self):
        for face in self.mcomplex.Faces:
            corner = face.Corners[0]
            face.InCenter = compute_incenter_of_triangle(
                [ corner.Tetrahedron.Class[v].IdealPoint
                  for v in simplex.ZeroSubsimplices
                  if simplex.is_subset(v, corner.Subsimplex) ])
        
    def _add_edge_centers(self):
        for edge in self.mcomplex.Edges:
            corner = edge.Corners[0]
            tet = corner.Tetrahedron
            info = self._edgeInfos[tet.Index][corner.Subsimplex]
            
            edge.MidPoint = compute_midpoint_of_triangle_edge_with_offset(
                [ tet.Class[info.perm.image(v)].IdealPoint
                  for v in [simplex.V0, simplex.V1, simplex.V2] ],
                info.offset)

    def _check_consistency(self, epsilon = 1e-6):
        for tet in self.mcomplex.Tetrahedra:
            for perm in t3m.Perm4.S4():
                v0 = perm.image(simplex.V0)
                v1 = perm.image(simplex.V1)
                v2 = perm.image(simplex.V2)
                v3 = perm.image(simplex.V3)

                edge = v0 | v1
                face = v0 | v1 | v2

                inRadius, inCenter = compute_inradius_and_incenter(
                    [ tet.Class[v].IdealPoint for v in [v0, v1, v2, v3]])
                inCenterTri = compute_incenter_of_triangle(
                    [ tet.Class[v].IdealPoint for v in [v0, v1, v2]])
                
                info = self._edgeInfos[tet.Index][edge]
                midpointEdge = compute_midpoint_of_triangle_edge_with_offset(
                    [ tet.Class[info.perm.image(v)].IdealPoint
                      for v in [simplex.V0, simplex.V1, simplex.V2] ],
                    info.offset)

                if not abs(inRadius - tet.InRadius) < epsilon:
                    raise Exception("InRadius %r %r" % (inRadius, tet.InRadius))

                if not inCenter.dist(tet.InCenter) < epsilon:
                    raise Exception("incenter %r %r %r" % (
                            inCenter,
                            tet.InCenter,
                            inCenter.dist(tet.InCenter)))

                if not inCenterTri.dist(tet.Class[face].InCenter) < epsilon:
                    raise Exception("incentertri %r %r %r" % (
                            inCenterTri,
                            tet.Class[face].InCenter,
                            inCenterTri.dist(tet.Class[face].InCenter)))
                
                if not midpointEdge.dist(tet.Class[edge].MidPoint) < epsilon:
                    raise Exception("midpoint %r %r %r" % (
                            midpointEdge,
                            tet.Class[edge].MidPoint,
                            midpointEdge.dist(tet.Class[edge].MidPoint)))


_transposition23 = t3m.Perm4([0,1,3,2])
            
def _find_perm_for_edge(e):
    for p in t3m.Perm4.A4():
        if e == p.image(simplex.E01):
            return p

def _invert_and_normalize(distances):
    geometric_mean = prod(distances) ** Rational((1, len(distances)))
    return [ geometric_mean / distance for distance in distances ]

def _doctest():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    _doctest()

from sage.all import RealIntervalField, ComplexIntervalField, prod, vector, matrix, arccosh, Infinity

from snappy.verify.upper_halfspace.finite_point import FinitePoint
from snappy.raytracing.hyperboloid_utilities import complex_and_height_to_R13_time_vector, PSL2C_to_O13

class PrecisionExperiment:
    def __init__(self, tiling_engine, bits_prec = 53):

        self.RIF = RealIntervalField(bits_prec)
        self.CIF = ComplexIntervalField(bits_prec)
        
        self.baseTetInCenter = FinitePoint(
            self.CIF(tiling_engine.baseTetInCenter.z),
            self.RIF(tiling_engine.baseTetInCenter.t))

        self.generator_matrices = {
            g : m.change_ring(self.CIF)
            for g, m
            in tiling_engine.mcomplex.GeneratorMatrices.items() }
    
        self.max_values = {}

        for tile in tiling_engine.all_tiles():
            if tile.word:
                matrix = prod(self.generator_matrices[g] for g in tile.word)
                tileCenter = FinitePoint(
                    self.CIF(tile.center.z),
                    self.RIF(tile.center.t))
                err = tileCenter.dist(self.baseTetInCenter.translate_PSL(matrix)).upper()
                l = len(tile.word)
                self.max_values[l] = max(err, self.max_values.get(l, err))

                #print("len=", l)
                #print(tile.center)
                #print(self.baseTetInCenter.translate_PSL(matrix))

def inner_prod(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b

    return a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3

def my_dist(a, b):
    RIF = a.base_ring()

    i = inner_prod(a,b)

    if not i > 0.9:
        raise Exception("Bad inner product", i)

    return arccosh(i.intersection(RIF(1, Infinity)))

class SO13PrecisionExperiment:
    def __init__(self, tiling_engine, bits_prec = 2 * 53):
        self.RIF = RealIntervalField(bits_prec)
        self.CIF = ComplexIntervalField(bits_prec)

        self.baseTetInCenter = vector(
            self.RIF,
            complex_and_height_to_R13_time_vector(
                tiling_engine.baseTetInCenter.z,
                tiling_engine.baseTetInCenter.t))

        self.generator_matrices = {
            g : matrix(self.RIF,
                       PSL2C_to_O13(m))
            for g, m
            in tiling_engine.mcomplex.GeneratorMatrices.items() }

        self.max_values = {}

        for tile in tiling_engine.all_tiles():
            if tile.word:
                m = prod(self.generator_matrices[g] for g in tile.word)
                tileCenter = vector(
                    self.RIF,
                    complex_and_height_to_R13_time_vector(
                        tile.center.z,
                        tile.center.t))

                #print("=====")
                #print(tileCenter)
                #print(m * self.baseTetInCenter)
                #print(inner_prod(m * self.baseTetInCenter, tileCenter).endpoints())

                err = my_dist(m * self.baseTetInCenter, tileCenter).upper()
                l = len(tile.word)
                self.max_values[l] = max(err, self.max_values.get(l, err))

if __name__ == '__main__':
    from snappy import Manifold
    from spineTilingEngine import get_tiling_engine
    import time

    t = get_tiling_engine(Manifold("m015"), 2.0, 1000)

    print(time.process_time())

    p = SO13PrecisionExperiment(t)

    for l in sorted(p.max_values.keys()):
        print(l, p.max_values[l])    

"""
SO13

# With cut_off 2.0, precision 53

1 3.14693332939488e-7
2 1.94607111421966e-6
3 8.14470808624817e-6
4 0.0000198918837387775
5 0.0000437288532913095
6 0.000177651761008844
7 0.000158143242844043
8 0.000264805666663187
9 0.000400416033450115
10 0.000424379709861500
11 0.000810140556132233
12 0.000843544180585911
13 0.000951147032856916
14 0.00169059297569973
15 0.00107394394002549
16 0.00272983953461852
"""

"""
# With cut_off 2.0, precision 106

1 3.270921281085800238169841633950e-15
2 1.958025535935146422497115696615e-14
3 8.202076444245559868708528727807e-14
4 2.056665304271714819705992280621e-13
5 4.395612855950487133117295347254e-13
6 1.837457195065786916202026100860e-12
7 1.642014914084237702968333803705e-12
8 2.775231119033088881471601358991e-12
9 4.219929997444452653216752991498e-12
10 4.405387882386405780612489783424e-12
11 8.354387324326215904772343210832e-12
12 8.686169163671624886197965679412e-12
13 9.801694308813672806298050836381e-12
14 1.714340993607170713978996529316e-11
15 1.105424563281960868139073605716e-11
16 2.812018996585752241579728815233e-11

"""

if __name__ == '__main__a':
    from snappy import Manifold
    from spineTilingEngine import get_tiling_engine
    import time

    t = get_tiling_engine(Manifold("m015"), 3.0, 1000)

    print(time.process_time())

    p = PrecisionExperiment(t)

    for l in sorted(p.max_values.keys()):
        print(l, p.max_values[l])

"""
PSL(2,C)

# With cut_off 3.0, precision 53

1 2.10734242554471e-8
2 2.10734242554471e-8
3 2.10734242554471e-8
4 2.10734242554471e-8
5 2.10734242554471e-8
6 2.10734242554471e-8
7 2.10734242554471e-8
8 2.10734242554471e-8
9 2.10734242554471e-8
10 2.10734242554471e-8
11 2.10734242554471e-8
12 2.10734242554471e-8
13 2.10734242554471e-8
14 2.10734242554471e-8
15 2.10734242554471e-8
16 2.10734242554471e-8
17 2.10734242554471e-8
18 2.10734242554471e-8
19 3.65002414998886e-8
20 4.21468485108941e-8
21 4.71216091538725e-8
22 4.71216091538725e-8
23 3.65002414998886e-8
24 5.16191365590357e-8
25 5.16191365590357e-8
26 5.57550398524693e-8
27 6.66400187462506e-8
28 5.57550398524693e-8
29 9.18569267238524e-8
30 8.68879540986613e-8
31 1.01064592348416e-7
32 8.94069671630860e-8
33 1.68587394043576e-7
34 3.44342589625457e-7
35 6.21219384243355e-7
36 5.20475184815312e-7
37 8.29393866662702e-7
38 8.92080637639480e-7
39 6.93184336996840e-7
40 1.24261747237373e-6
41 7.69396462644653e-7

"""

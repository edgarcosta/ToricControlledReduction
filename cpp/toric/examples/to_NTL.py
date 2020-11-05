from sage.all import Integer, GF, LaurentPolynomialRing, Matrix, Polyhedron, PolynomialRing, ZZ
from sage.all import vector

# makes the vector v immutable and returns it
def immutable(v):
    v = vector(v);
    v.set_immutable()
    return v

def vector_to_NTL(v):
    return str(list(v)).replace(", ", " ")
def matrix_to_NTL(m):
    return "[ " + Matrix(m).str().replace("\n","") + " ]"
def map_to_NTL(m):
    return matrix_to_NTL(Matrix(map(vector, m.keys()))) + "\n" + vector_to_NTL(m.values())


def input_to_NTL_Fp(f, p):
    #p = vector(f.values()).base_ring().characteristic()
    n = len(f.keys()[0])
    P = Polyhedron(f.keys() + [[0]*n])
    AP = Matrix([l[1:] for l in P.inequalities_list()])
    bP = vector(l[0] for l in P.inequalities_list())
    output = "%s \n%s \n%s \n%s \n" % (p, map_to_NTL(f), matrix_to_NTL(AP), vector_to_NTL(bP))
    return output


def ZZ_pE_to_NTL(ola):
    return vector_to_NTL(vector(ola.polynomial().list()))

def vector_ZZ_pE_to_NTL(ola):
    res = "[ ";
    for v in ola:
        res += ZZ_pE_to_NTL(v) + " "
    res += " ]"
    return res
def matrix_ZZ_pE_to_NTL(ola):
    res = "[ ";
    for v in ola.rows():
        res += vector_ZZ_pE_to_NTL(v) + " "
    res += " ]"
    return res
def map_ZZ_pE_to_NTL(m):
    return matrix_to_NTL(Matrix(map(vector, m.keys()))) + "\n" + vector_ZZ_pE_to_NTL(m.values())

def input_to_NTL_Fq(f):
    Fq = f.values()[0].parent()
    assert Fq.degree() > 1;
    p = Fq.characteristic()
    fE = vector(Fq.polynomial().list())
    n = len(f.keys()[0])
    P = Polyhedron(f.keys() + [[0]*n])
    AP = Matrix([l[1:] for l in P.inequalities_list()])
    bP = vector(l[0] for l in P.inequalities_list())
    ffrob = {}
    for v, fv in f.iteritems():
        ffrob[v] = fv.frobenius();
        
    output = "%s \n%s \n%s \n%s \n%s \n%s \n" % (p, vector_to_NTL(fE), map_ZZ_pE_to_NTL(f), map_ZZ_pE_to_NTL(ffrob), matrix_to_NTL(AP), vector_to_NTL(bP))
    return output

def matrix_ZZq_to_NTL(F):
    Zq = F.base_ring();
    ZZX = PolynomialRing(ZZ, "X")
    p = Zq.prime()

    Fx = [[ ZZX(sum( [ZZX(ci)*p^i for i,ci in enumerate(elt.list())])) for elt in row] for row in F.rows()]
    res = "[ "
    for row in Fx:
        res += "[ "
        for elt in row:
            res += vector_to_NTL(vector(ZZX(elt).list())) + " ";
        res += "] "
    res += "] "
    return res

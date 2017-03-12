sys.path.append("../../../sage")
from toric import hypersurface
from sage.all import Integer, GF, LaurentPolynomialRing, Matrix, Polyhedron, PolynomialRing, ZZ
from sage.all import vector

def vector_to_NTL(v):
    return str(list(v)).replace(", ", " ")
def matrix_to_NTL(m):
    return "[ " + Matrix(m).str().replace("\n","") + " ]"
def map_to_NTL(m):
    return matrix_to_NTL(Matrix(map(vector, m.keys()))) + "\n" + vector_to_NTL(m.values())


def input_to_NTL_Fp(f):
    p = vector(f.values()).base_ring().characteristic()
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

# CURVES
R2 = LaurentPolynomialRing(ZZ, ("x","y"))
x, y = R2.gens();

data = [];
# 11.a elliptic curve
# http://www.lmfdb.org/EllipticCurve/Q/11/a/1
data += [[17, y**2 + y - (x ** 3 - x**2 - 7820*x - 263580), [17, 2, 1],  "elliptic curve 11.a p = 17", "ecp"]]
data += [[17**2, y**2 + y - (x ** 3 - x**2 - 7820*x - 263580), [289, 30, 1],  "elliptic curve 11.a p = 17","ecq"]]

# genus 2 curve
# http://www.lmfdb.org/Genus2Curve/Q/11664/a/11664/1
data += [[31, (y - 16)**2 + (y - 16) + x**6, [961, 0, 46, 0, 1], "genus 2 curve 11664.a  p = 31","g2p"]]
data += [[31**2, (y - 16)**2 + (y - 16) + x**6, [923521, 88412, 4038, 92, 1], "genus 2 curve 11664.a q = 31**2","g2q"]]

# C_{3,4} curve
data += [[43, -y**3 +  x**4 + x + 1, [79507, 27735, 6579, 1258, 153, 15, 1], "C_{3,4} curve p = 43","c34"]]

# dense quartic curve on P^2
data += [[17, -2*x**4 - 3*x**3*y - x**2*y**2 + x*y**3 + y**4 - x*y**2 + 7*y**3 + 5*x**2 - 2*x*y - 2*y**2 + 6*y + 1, [4913, 867, 221, 27, 13, 3, 1], "plane curve, genus 3 p = 17","pg3"]]

#genus 9 curve
data += [[37, y**4 - x**3 + 1/(x*y), [37**9, 0, 9*37**8,0, 36*37**7,0,84*37**6,0,126*37**5, 0, 126*37**4,0, 84*37**3, 0,36*37**2,0, 9*37, 0, 1], "genus 9 curve p = 37","g9"]]

# hyperelliptic curve of genus 10
data += [[13, -y**2 +  x**22 + x + 1, [137858491849, 31813498119, 8973037931, 3074677333, 323396203, 24876631, 10995985, 12881011, 3433911, 999635, 664001, 76895, 20319, 5863, 385, 67, 67, 49, 11, 3, 1], "hyperelliptic curve of genus 10 p = 13","g10"]]
data += [[103, -y**2 +  x**22 + x + 1,[134391637934412192049, 5219092735316978332, -152012409766513932, -105769152426538820, 186272158258524, -534123936992982, 124637044805304, 8798421444054, 969104502946, -108512162008, -4925241640, -1053516136, 91347394, 8051802, 1107384, -46074, 156, -860, -12, 4, 1], "hyperelliptic curve of genus 10 p = 103","g10_103"]]

def curves_Fp():
    print """
// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    bool val = true;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
"""
    for q, f, expected, desc, short in data:
        Fq = GF(Integer(q))
        fbar = R2.change_ring(Fq)(f).dict()

        if Fq.degree() == 1:
            print "// %s"  % desc
            print "const char %s[] = %s;" % (short, repr(input_to_NTL_Fp(fbar)).replace("'","\""))
            print "const char %s_zeta[] = %s;" % (short, repr(vector_to_NTL(expected)).replace("'","\""))
            print "cout << \"Testing: %s\" <<endl;" % desc
            print "user_time = get_cpu_time();"
            print "get_timestamp(&wtime1);"
            print "if( test_Fp(%s, %s_zeta, 0))\ncout <<\"PASS\"<<endl;\nelse\nval = false;" % (short, short)
            print "get_timestamp(&wtime2);"
            print "wall_time = timestamp_diff_in_seconds(wtime1,wtime2);"
            print "user_time = get_cpu_time() - user_time;"
            print "printf(\"Time: CPU %.2f s, Wall: %.2f s\\n\\n\\n\", user_time, wall_time );"
            print 

    print "return not val;\n}\n"

def curves_Fq():
    print """
// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    bool val = true;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
"""
    
    for q, f, expected, desc, short in data:
        Fq = GF(Integer(q))
        fbar = R2.change_ring(Fq)(f).dict()
        xq, yq = R2.change_ring(Fq)(f).parent().gens()
        fbar = R2.change_ring(Fq)(f).subs(xq = xq * Fq.gen(), yq = (1 + Fq.gen())*yq).dict()
        if Fq.degree() != 1:
            H = hypersurface(Fq, fbar, 0)
            H.zeta_function()
            F = Matrix(H.padic.Zq, H.frob_matrix)
            print "// %s"  % desc
            print "const char %s[] = %s;" % (short, repr(input_to_NTL_Fq(H.fbar)).replace("'","\""))
            print "const char %s_F[] = %s;" % (short, repr(matrix_ZZq_to_NTL(F)).replace("'","\""))
            print "user_time = get_cpu_time();"
            print "get_timestamp(&wtime1);"
            print "cout << \"Testing: %s\" <<endl;" % desc
            print "if( test_Fq(%s, %s_F, 0))\ncout <<\"PASS\"<<endl;\nelse\nval = false;" % (short, short)
            print "get_timestamp(&wtime2);"
            print "wall_time = timestamp_diff_in_seconds(wtime1,wtime2);"
            print "user_time = get_cpu_time() - user_time;"
            print "printf(\"Time: CPU %.2f s, Wall: %.2f s\\n\\n\\n\", user_time, wall_time );"
            print 
    print "return not val;\n}"


R3 = LaurentPolynomialRing(ZZ, ("x","y", "z"))
x, y, z = R3.gens()
# dwork K3 surface
data2 = []
data2 += [[43, y**4 +  x**4 + z**4 + x*y*z + 1, [-20083415214428110320965436874242043, 119479484780264582763991241544977, 100534534762252228203104681046302, -834165241103179642680896617306, -225003300376177017932433136047, 2487646096419799283278902381, 296132512235307123091919592, -4247548915225822821120696, -253423978516239259535142, 4644610109877171089586, 147044525659586030196, -3419640131618279772, -58417624987449798, 1723873631640594, 15626409457128, -589209839544, -2676921183, 130948029, 262558, -17114, -11, 1], "dwork K3 surface p = 43", "dwork"]]


#small volume polyhedron K3 surfaces
data2 += [[43, y**1*z**3 + y**2*z + y**3 + x**4 + z**4 + 1, [-20083415214428110320965436874242043, -467056167777397914441056671494001, 18944950018012354560886560498675, 693179566550529562227787329874, 1229526231563808841160836809, -346300669871977167148659463, -12191083961782853121299085, -96222790959808283973608, 6233702404372961918556, 157050641970636637483, 1015739668724359703, -23621852761031621, -1975305846914569, -42403703332692, 353996330744, 24256393095, 372649309, -715563, -218182, -3225, 43, 1], "small volume polyhedron K3 surface #1 p = 43", "smallvol1"]];


data2 += [[43, y**2*x**2 + y**3 + x**4 + z**4 + 1, [-20083415214428110320965436874242043, 467056167777397914441056671494001, 10609172010086918554096473879258, -246724930467137640792941253006, 12158648289908776318146052889, -282759262556018053910373323, -6501911446284188331359512, 151207242936841589101384, -1838143016674078514446, 42747512015676244522, 1988256372822150908, -46238520298189556, -537657212769646, 12503656110922, -556279948312, 12936742984, 304273289, -7076123, 77658, -1806, -43, 1],  "small volume polyhedron K3 surface #2 p = 43","smallvol2"]];

#degree 2 K3 surface
# computed using Magma
# f6 := x^6 + y^6 + z^6 + x^2*y^2*z^2; wf1, wf2 := WeilPolynomialOfDegree2K3Surface(PolynomialRing(GF(43),3)!f6); wf1*wf2;
data2 += [[43, -y**2 + x**6  + z**6 + 1 + x**2*z**2, [20083415214428110320965436874242043, -5441747443173868724069055637639407, 635034724603774124880917507915586, -40122173216441668728947351381690, 1326658803857349739612542916911, -8492309087777935484296942611, -1090547874399484315578027240, 41939390745482153486211144, -406389445208160401998170, -10461988918967025235406, 205654865693212913484, 4782671295190997988, -131585758725232058, -2764395274958190, 154291829299992, -2169844618680, -9138455073, 772092477, -12628670, 108102, -501, 1], "degree 2 K3 surface p = 43","degree2"]];

# JV HGM
# (x^3*y + y^4 + z^4 - 12*x*y*z*w + w^4).subs(w = x, x = 1, y = y + 1)
data2 += [[71, x^3*y + y^4 + z^4 - 12*x*y*z + 1, [-2102085018129621311776010144838961, 1667990492465480112498321876484, 1180514397836484164127333440751, -992654528346843947132506572, -248163632086710986783126643, 131277461925126489073928, 54775983232843270969227, 13020974204039524804, -18981279297437898834, 0, 3765379745573874, -512401135684, -427602356187, -203293448, 76235043, 60492, -14271, -4, 1],"hypergeometric C2F2, small volume, p = 71 ","c2f20"]]
data2 += [[71, x^4 + y^4 + z^4 + 4*y^3 - 12*x*y*z + 6*y^2 - 12*x*z + 5*y + 2,[752359350923790893319063566949457370471, 9999618431242608580118480258998937477, -580175465003775482050178799975553236, -7579348797763141631192374606743296, 179308977940039112093569601722030, 2314043147997884381423728237094, -38878969026918054596023174070, -519629013665131925376875506, 10748342195119731292927855, 151385101339714525252505, -2682320686032142109624, -37779164592002001544, 422968562965973455, 5957303703746105, -57132726628766, -847987794970, 10012202314, 153901730, -1290496, -19596, 67, 1],"hypergeometric C2F2, sparse, p = 71 ","c2f22"]]
# ??
data2 += [[71, 35*x^4 + 68*x^3*y + 48*x^2*y^2 + 56*x*y^3 + 40*y^4 + 25*x^3*z + 16*x^2*y*z + 48*x*y^2*z + 43*y^3*z + 38*x^2*z^2 + 61*x*y*z^2 + 4*y^2*z^2 + 69*x*z^3 + 67*y*z^3 + 57*z^4 + 12*x^3 + 63*x^2*y + 68*x*y^2 + 59*y^3 + 28*x^2*z + 11*x*y*z + 63*y^2*z + 5*x*z^2 + 64*y*z^2 + 40*z^3 + 28*x^2 + 29*x*y + 16*y^2 + 43*x*z + 64*y*z + 12*z^2 + 63*x + 56*y + 29*z + 11,[752359350923790893319063566949457370471, 9999618431242608580118480258998937477, -580175465003775482050178799975553236, -7579348797763141631192374606743296, 179308977940039112093569601722030, 2314043147997884381423728237094, -38878969026918054596023174070, -519629013665131925376875506, 10748342195119731292927855, 151385101339714525252505, -2682320686032142109624, -37779164592002001544, 422968562965973455, 5957303703746105, -57132726628766, -847987794970, 10012202314, 153901730, -1290496, -19596, 67, 1],"hypergeometric C2F2, random change of variables, p = 71","c2f2"]]

# Sturmfels example
Sturmfels2 = 22*x^4 + 44*x^3*y + 75*x^2*y^2 + 57*x*y^3 + 21*y^4 - 16*x^3*z + 43*x^2*y*z + 48*x*y^2*z + 49*y^3*z - 27*x^2*z^2 + 5*x*y*z^2 + 64*y^2*z^2 - 40*x*z^3 - 22*y*z^3 - 31*z^4 - 56*x^3 + 27*x^2*y + 51*x*y^2 + 48*y^3 - 116*x^2*z - 64*x*y*z + 51*y^2*z - 125*x*z^2 - 25*y*z^2 - 87*z^3 - 57*x^2 - 32*x*y + 38*y^2 - 126*x*z - 49*y*z - 157*z^2 - 72*x + 21*y - 108*z - 49;
data2 += [[107, Sturmfels2, [-4140562374860211619063098135818149338786107, -43759983173909127950618820371560491745403, 804423881624343280470169093105074305662, 14435785514109985045742059041428530499, -12694304118831251292805753466971243, -1583684104370601660958596498596282, -9850000453768762162333255111918, 18555806501448223479465201439, 774753906602523228175267544, 10061192331993210798407354, 51736080696715567981487, -483514772866500635341, -8212929939596578078, -55238882790857992, -115556237092373, 5357747611274, 75239690974, 52676849, -5232193, -25466, 121, 1],"dense, rank 2, p = 107","r2_107"]]
data2 += [[223, Sturmfels2, [20625387851152557293083391544660009040933060233823, 14516450658375183600271847494683189214194074045, -684439954168106030096097357853095916825550192, -483738513133044328153315055998727373463174, 4450663604672357743176472217709027171959, -11740078091987226966965107406249082490, 87994179164798562148297089463878453, 529496229691210542228455906091331, -2646649344376391559509290896175, -10444176784983069831247425958, 52917324575943572029498926, 237297419623065345423762, -941802036543272594074, -4799235816058258225, 19307623910855773, 64522350049131, -173108140870, 1319658473, -2884282, -82064, 35, 1], "dense, rank 2, p = 211", "r2_211"]]
data2 += [[1009, Sturmfels2, [-1207020068576253752077404803423021173724331179038253595707180209, 735061790287096337411258021829572625075102404429232280475180, 598079304415313845938697421770600828292507950788834730589, 758107537417407017989551836171239082478562633671283419, -298922374803585575685640183918239564901977018082819, 93795358123965458388872713460307595815123791218, -289078976064337932862318909477851010815545055, -785348352159864027947534610952478114848971, -161458902868717003410496544813572449745, 116979212631721420057500657873884554, 519523589581200090948564920340475, -514889583331219118878657007275, -113876784378716241408664426, 154385378152419311145505, 737605545460107983259, 266683417330332495, -84992091650002, 266056125811, -662770731, -513581, -620, 1],"dense, rank 2, p = 1009","r2_1009"]]

def k3_surfaces():
    print """
// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    bool val = true;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;
"""
    for q, f, expected, desc, short in data2:
        Fq = GF(Integer(q))
        fbar = R3.change_ring(Fq)(f).dict()

        if Fq.degree() == 1:
            print "// %s"  % desc
            print "const char %s[] = %s;" % (short, repr(input_to_NTL_Fp(fbar)).replace("'","\""))
            print "const char %s_zeta[] = %s;" % (short, repr(vector_to_NTL(expected)).replace("'","\""))
            print "user_time = get_cpu_time();"
            print "get_timestamp(&wtime1);"
            print "cout << \"Testing: %s\" <<endl;" % desc
            print "if( test_Fp(%s, %s_zeta, 0))\ncout <<\"PASS\"<<endl;\nelse\nval = false;" % (short, short)
            print "get_timestamp(&wtime2);"
            print "wall_time = timestamp_diff_in_seconds(wtime1,wtime2);"
            print "user_time = get_cpu_time() - user_time;"
            print "printf(\"Time: CPU %.2f s, Wall: %.2f s\\n\\n\\n\", user_time, wall_time );"
            print
    print "return not val;\n}"

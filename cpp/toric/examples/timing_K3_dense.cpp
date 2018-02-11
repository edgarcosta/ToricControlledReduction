
// Copyright 2018 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#define verbose 0

int main()
{

    int examples_length = 11;

    char examples[11][3][buffer_length] = {
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 127 = 2^7 - 1",
        "127 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[-151313661355466579537756144585602921111718527 -1275879344307982814627989067123938078690168 9676945003725437049343826296505088484739 102370759190700061350274293934296782672 73278997273228390372422543975874576 -288499989264678702253632062897144 -8802657152760865915219088533279 -129681004536482149627306866886 739424679681272099422602282 -3604244782848345441504998 -90597309805091663126867 713364644134580024621 1759556090266490906 -22380786126951126 243360629051962 1024189451617 2081157128 -32774128 -2838704 -16637 136 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 239 = 2^8 - 17",
        "239 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[-88380391271935664108453035995755078242378693829039 -126874390929758310549415257990089746605890178638 -97107657105131883775958203380515675022367215 -541744251632534916462807271299947977809580 18247034416911740910148947840855988373921 -26555625856884469216152734714725833544 33730248025032574467617789162458097 1311685746481699917766086804742402 764183211233446273789892211862 1162697925041378887470357112 32229597294557050751008853 -134851871525343308581627 -85167361822274135048 -979959184147505738 -29447232996444638 -13256791508383 182717283896 -2197958959 1142420 3585 82 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 509 = 2^9 - 3",
        "509 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[693543058066147973304560948699228887716627117992687269709 -294462876039834171797629715636867148300450372583074790 1036063343094254899877354062876635617345753230913837 6106463711753958152518786225206103835043732205779 -5216961607879712535436019865360793984443267137 -2233339576657242796850596748489863110920344 56178237316279928618505398142970525208037 -102519245088113900858260091472650916114 -116130049835205749746520219017144704 231718235488960392114557165892510 -853287139904438298502550442131 -1676399096079446558944106959 1757142024868332149861190 -3399035502238877364096 -11581929664431341706 24496727285430633 -3758886015416 -33891162853 153116871 100273 -110 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 1021 = 2^10 - 3",
        "1021 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[-1547179077535659880232325182481393622240140032528876996515545421 0 236946862252738349042554301324957357865095123055944792443 52679220440391766249704122320973797137466092746014417 1084904399571823970730068749919729310769744028281738 -285451557788386316508318082142895626233261535909 -16052461436440659069117044391766902299552172 690470713108409021337306757628236819446947 7699458020377488543292639291704231846 437383511441620309021521135082120908 615499104216089011535795894100500 -602839475236130275745147790500 -410946400356852763877201628 -6939552169219709613246 -596987266271431107767 13314043037869212 227116796842529 -828050499058 -38570317 -166423 0 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 2039 = 2^11 - 9",
        "2039 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[-3146021898873573193142129399165443470209827369368071531647907531226839 -231552095841563613774047466301342964204920955306450620137399114654 37482749265143308412714796911909976671726523135382127971672181 -350185295578475788393957770303414094058929248527356227899996 71678871249082613078327114396513699473907500866030046883 -26091836378937760613664281370267434398133574937629004 -47599989180859988019273258537283789971941548887897 2516644293787060559196327630689782396912446841 -1853963407014174929361717675654510847042439 -1691865934452915846854891465089819885252 561449491188822425853352839847963652 -275355316914576962164469269175068 199578736430465882900547153308 52603546985384737118183761 -17175178352907922490879 78136116457934019183 10301848595802836 -6807179811157 7999070404 -205939 306 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 4093 = 2^12 - 3",
        "4093 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[7126505711707377646777607190509780940799469487793902854463253107062275702093 935445257645195302860100406325650844163095245633021522564263931180671443 -1021344747232498596046445682777889161753969746241595334822696282744123 -56321109391388238267222472955135278400131333818620861963439289474 79683465053267893467877863417182331770476067539530790881440044 -2713183708559347881079396698755575159333722060196052867030 -3370708965614707275616315736237675865765578538082387998 656326941044605650372434459320364820769866660225822 15385449013733883131825745640727025363625834296 -54569820928767911560516600375292363565098328 7899992965911265727052213936866395996763 1930122884415163871744982637885755191 -795842805765240591528052802947704 13393727293625103025856398872 34105782576629463167267046 -10455518008766585413686 -502365734845789790 880695000233308 -37157375482 -40221911 2199 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 8191 = 2^13 - 1",
        "8191 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[15138242034681741844672511532539354571077308957635664919842369432269759886351572991 -162681009024029620696377002587829098848598843117089594692531169933431093091631 -159741506370756674415373566029021898332482691963084979881675255433929447079 -5854994776725159361988966655631357407528393901257817721699490614974771 73903212987688094164096273695113059489517669589330663066186547380 152780825286882154541003976059567143697045440733243709915853848 6651912183192907377138186444020944071140128275874612102207 -1731035884456616495613993502954644995717097581264049907 88747463943171616364119744170871963002896688907213 37704495924396389365375133959674295557366461126 -7108662591030058129708819747313778815390829 -867862604203401065768382340045632867219 68609200925361479490039679291886906 2406973571701293156593961229363 -699757544236680226176125197 40078750527386074700737 13720270268529380328 98919812136780 -116808009421 -47499609 -721 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 16381 = 2^14 - 3",
        "16381 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[-31706522829727082208656711649602609814876533431553847447248808771529041219123982028827981 -1549304337103654092546305217671992213902349165699543373078790159117353069354242934552 34457409726681872133209302742566653480054115789972893291554370914284110366935857 6134800456531653515420589061591311882195653784525792591673179007184263147452 223086006769098106070676146349753891565085035360425464723306274806962219 -15013425124813629546618578987677433361783835442037304383911142696449 -1184787166457638275587506468280877942567706378282168426832746267 -13080854613026545841147261537958707722515727860704840473599 -204581087238880216759541838066632397153654561131974628 -121835392835282913977305765303159303214405562859626 -987785056897111769852980345355568464608278710 60300656669135691951222779156069132812910 27717383894016634624899250107548644386 173445657376879852826090014029428 41328842914542039446069857059 13950078401310058912671127 658772207725916737229 -36479341953995959 -3738473327052 -78252037 13112 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 32749 = 2^15 - 19",
        "32749 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[65941513015777616542367829443641394522344565552758587358833918487147237226091144462310387562749 -1254890483980030405352832682601338912462891430436552584815852810317252450559085767061977090 -80985011569764455319137609082948769032356142038097070901232700607417257808137139371136 390919581597628106175948650288812457702258914251187131687613443676442275244971431 81455499467481183763692723031034166563265893604711264273036415298825507278532 -352895538960649384611338034948176368618011515601606997668124903286350898 -21629876131677640967715637860256950372717262564987015164833739721252 -253882930466748958985877839815121405122735146974125884765515406 -5582205826393660740628925394568825156286477684678693979668 85087706925468840775175653614532760889121362706178419 7255314202646153030124468335121217928318555810113 221543076205262848640400266729402971947801637 2422549669250780633123172318795842029831 -148188884740075667838349845684135332 -6284159507546339556819437159094 -499196628283707933687501748 -7593949292850668394602 1634352941824472468 7313357049819 -1412660864 -20410 1]"
      },
      {
        "K3_dense, p = 49999",
        "49999 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[-476636926645927885787432991384108929249530530805399402199082512717902131359025093916249475001049999 -2500728000074907856279476482545100297408593201582073426001477183409947266488669892922460186884 184256012408475483793860116543740910302325922863751886002900564301595946480481974012948319 -2102327682448277045914168604168547781429195530138513421428408722968015627923430222435 -10802816924902762912261069131937077622801549261125583150557254262124594334407082 1150930300424581580773238828239784940630155291261607848239349472215793712275 -18505559853424615308311885088785246973027651234738170082899598770430328 -115387475524589453632971942948793366615140560806764406281145740545 7558488461297436653444316178328398263303304662609536419830967 -32278054225855822798358520579593862886371590853636593388 -1095874550978635286827866124731331641262694388011224 21917929378160268941936161217850989845050788776 258239927892634292892703716462837355206612 -24189581937589163604228999243338419033 147716647761382612854559663509455 9476552287909241807581969672 -235762388065867455037725 885196886062292918 68909743527565 -2415901681 13116 1]"
      },
      // f = -9*x^4 - 10*x^3*y - 9*x^2*y^2 + 2*x*y^3 - 7*y^4 + 6*x^3*z + 9*x^2*y*z - 2*x*y^2*z + 3*y^3*z + 8*x^2*z^2 + 6*y^2*z^2 + 2*x*z^3 + 7*y*z^3 + 9*z^4 + 8*x^3 + x^2*y - 8*x*y^2 - 7*y^3 + 9*x^2*z - 9*x*y*z + 3*y^2*z - x*z^2 - 3*y*z^2 + z^3 - x^2 - 4*x*y - 3*x*z + 8*y*z - 6*z^2 + 4*x + 3*y + 4*z - 5
      {
        "K3_dense, p = 65521 = 2^16 - 15",
        "65521 \n[ [0 1 1][1 0 3][1 0 0][1 0 1][1 3 0][0 2 1][2 2 0][4 0 0][0 2 2][0 1 2][2 0 1][0 0 4][1 2 0][2 0 0][0 0 3][1 2 1][2 0 2][0 0 2][0 3 0][3 1 0][0 0 1][0 3 1][0 0 0][1 1 1][0 1 3][1 0 2][1 1 0][2 1 0][2 1 1][0 4 0][3 0 0][3 0 1][0 1 0] ]\n[8 2 4 -3 2 3 -9 -9 6 -3 9 9 -8 -1 1 -2 8 -6 -7 -10 4 3 -5 -9 7 -1 -4 1 9 -7 8 6 3] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n",
        "[139312748238933579349197720901032296262487215432324291509151053870310249151349047325114401620183199921 -690267988146793675007477063943957332997169510957303683954470880149739346360108893718982666622151 16209963778320057035488218061438124511269364339359992490754670277885936463095019599941622969 -179633914886786405367723761177061380541619482227676526303412274846696773944897364319124 -3063848506776027265309930985723725779346757277594543300324021376988008432778756797 -45870360763660865374665616097530120878949726291230705566301072406111892580851 -571280818772745446474776437985613119005313374268504963487386128028186698 -1940437552283966882584932571541069366546272095834468963488971876491 94787146073838930053879692323211266877612803533883655699770022 1349217282544396847575776969650188096941025982209767845962 30271400952455122870951473883464895092252087779252760 462010667609699529478357685069899651901712241560 4796675557302074852359008659593380218944042 78495879844830354927455632569675634822 -374313982781847842078474210351451 -25669947966961637861114414058 -480116318613123820406531 -7469999366120364877 -102018886243924 2144436809 -21271 1]"
      },

    };



    return not run_examples(examples, examples_length, verbose);
}

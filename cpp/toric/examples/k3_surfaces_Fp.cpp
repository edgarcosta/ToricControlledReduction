// Copyright 2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"

int main()
{
    bool val = true;
    timestamp_type wtime1, wtime2;
    double wall_time, user_time;

    // dwork K3 surface p = 43
    const char dwork[] = "43 \n[ [1 1 1][0 0 0][0 0 4][0 4 0][4 0 0] ]\n[1 1 1 1 1] \n[ [ 1  0  0][ 0  1  0][-1 -1 -1][ 0  0  1] ] \n[0 0 4 0] \n";
    const char dwork_zeta[] = "[-20083415214428110320965436874242043 119479484780264582763991241544977 100534534762252228203104681046302 -834165241103179642680896617306 -225003300376177017932433136047 2487646096419799283278902381 296132512235307123091919592 -4247548915225822821120696 -253423978516239259535142 4644610109877171089586 147044525659586030196 -3419640131618279772 -58417624987449798 1723873631640594 15626409457128 -589209839544 -2676921183 130948029 262558 -17114 -11 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: dwork K3 surface p = 43" <<endl;
    if( test_Fp(dwork, dwork_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // small volume polyhedron K3 surface #1 p = 43
    const char smallvol1[] = "43 \n[ [0 3 0][0 0 0][0 1 3][0 2 1][0 0 4][4 0 0] ]\n[1 1 1 1 1 1] \n[ [  1   0   0][  0   1   0][ -1  -1  -1][ -9 -12  -8][  0   0   1] ] \n[0 0 4 36 0] \n";
    const char smallvol1_zeta[] = "[-20083415214428110320965436874242043 -467056167777397914441056671494001 18944950018012354560886560498675 693179566550529562227787329874 1229526231563808841160836809 -346300669871977167148659463 -12191083961782853121299085 -96222790959808283973608 6233702404372961918556 157050641970636637483 1015739668724359703 -23621852761031621 -1975305846914569 -42403703332692 353996330744 24256393095 372649309 -715563 -218182 -3225 43 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: small volume polyhedron K3 surface #1 p = 43" <<endl;
    if( test_Fp(smallvol1, smallvol1_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // small volume polyhedron K3 surface #2 p = 43
    const char smallvol2[] = "43 \n[ [0 3 0][0 0 0][0 0 4][2 2 0][4 0 0] ]\n[1 1 1 1 1] \n[ [ 0  0  1][-2 -4 -3][ 1  0  0][ 0  1  0][-1 -1 -1] ] \n[0 12 0 0 4] \n";
    const char smallvol2_zeta[] = "[-20083415214428110320965436874242043 467056167777397914441056671494001 10609172010086918554096473879258 -246724930467137640792941253006 12158648289908776318146052889 -282759262556018053910373323 -6501911446284188331359512 151207242936841589101384 -1838143016674078514446 42747512015676244522 1988256372822150908 -46238520298189556 -537657212769646 12503656110922 -556279948312 12936742984 304273289 -7076123 77658 -1806 -43 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: small volume polyhedron K3 surface #2 p = 43" <<endl;
    if( test_Fp(smallvol2, smallvol2_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // degree 2 K3 surface p = 43
    const char degree2[] = "43 \n[ [6 0 0][0 0 0][0 0 6][0 2 0][2 0 2] ]\n[1 1 1 42 1] \n[ [ 0  1  0][ 1  0  0][ 0  0  1][-1 -3 -1] ] \n[0 0 0 6] \n";
    const char degree2_zeta[] = "[20083415214428110320965436874242043 -5441747443173868724069055637639407 635034724603774124880917507915586 -40122173216441668728947351381690 1326658803857349739612542916911 -8492309087777935484296942611 -1090547874399484315578027240 41939390745482153486211144 -406389445208160401998170 -10461988918967025235406 205654865693212913484 4782671295190997988 -131585758725232058 -2764395274958190 154291829299992 -2169844618680 -9138455073 772092477 -12628670 108102 -501 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: degree 2 K3 surface p = 43" <<endl;
    if( test_Fp(degree2, degree2_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // hypergeometric C2F2, small volume, p = 71 
    const char c2f20[] = "71 \n[ [1 1 1][3 1 0][0 0 0][0 4 0][0 0 4] ]\n[59 1 1 1 1] \n[ [ 1  0  0][-1  3  0][-1 -1 -1][ 0  0  1] ] \n[0 0 4 0] \n";
    const char c2f20_zeta[] = "[-2102085018129621311776010144838961 1667990492465480112498321876484 1180514397836484164127333440751 -992654528346843947132506572 -248163632086710986783126643 131277461925126489073928 54775983232843270969227 13020974204039524804 -18981279297437898834 0 3765379745573874 -512401135684 -427602356187 -203293448 76235043 60492 -14271 -4 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: hypergeometric C2F2, small volume, p = 71 " <<endl;
    if( test_Fp(c2f20, c2f20_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // hypergeometric C2F2, sparse, p = 71 
    const char c2f22[] = "71 \n[ [0 3 0][0 0 0][0 2 0][1 1 1][0 1 0][0 0 4][4 0 0][1 0 1][0 4 0] ]\n[4 2 6 59 5 1 1 59 1] \n[ [ 0  1  0][ 1  0  0][ 0  0  1][-1 -1 -1] ] \n[0 0 0 4] \n";
    const char c2f22_zeta[] = "[752359350923790893319063566949457370471 9999618431242608580118480258998937477 -580175465003775482050178799975553236 -7579348797763141631192374606743296 179308977940039112093569601722030 2314043147997884381423728237094 -38878969026918054596023174070 -519629013665131925376875506 10748342195119731292927855 151385101339714525252505 -2682320686032142109624 -37779164592002001544 422968562965973455 5957303703746105 -57132726628766 -847987794970 10012202314 153901730 -1290496 -19596 67 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: hypergeometric C2F2, sparse, p = 71 " <<endl;
    if( test_Fp(c2f22, c2f22_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // hypergeometric C2F2, random change of variables, p = 71
    const char c2f2[] = "71 \n[ [0 3 0][3 1 0][0 0 0][2 0 0][0 2 2][0 2 1][0 1 0][4 0 0][1 0 2][0 1 1][1 1 2][1 1 0][1 3 0][2 2 0][2 0 1][1 0 1][1 0 0][2 0 2][0 2 0][0 1 2][1 2 1][3 0 0][0 0 3][0 3 1][0 4 0][1 0 3][0 0 1][1 1 1][0 1 3][0 0 2][3 0 1][1 2 0][2 1 0][2 1 1][0 0 4] ]\n[59 68 11 28 4 63 56 35 5 64 61 29 56 48 28 43 63 38 16 64 48 12 40 43 40 69 29 11 67 12 25 68 63 16 57] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n";
    const char c2f2_zeta[] = "[752359350923790893319063566949457370471 9999618431242608580118480258998937477 -580175465003775482050178799975553236 -7579348797763141631192374606743296 179308977940039112093569601722030 2314043147997884381423728237094 -38878969026918054596023174070 -519629013665131925376875506 10748342195119731292927855 151385101339714525252505 -2682320686032142109624 -37779164592002001544 422968562965973455 5957303703746105 -57132726628766 -847987794970 10012202314 153901730 -1290496 -19596 67 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: hypergeometric C2F2, random change of variables, p = 71" <<endl;
    if( test_Fp(c2f2, c2f2_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // dense, rank 2, p = 107
    const char r2_107[] = "107 \n[ [0 3 0][3 1 0][0 0 0][2 0 0][0 2 2][0 2 1][0 1 0][4 0 0][1 0 2][0 1 1][1 1 2][1 1 0][1 3 0][2 2 0][2 0 1][1 0 1][1 0 0][2 0 2][0 2 0][0 1 2][1 2 1][3 0 0][0 0 3][0 3 1][0 4 0][1 0 3][0 0 1][1 1 1][0 1 3][0 0 2][3 0 1][1 2 0][2 1 0][2 1 1][0 0 4] ]\n[48 44 58 50 64 51 21 22 89 58 5 75 57 75 98 88 35 80 38 82 48 51 20 49 21 67 106 43 85 57 91 51 27 43 76] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n";
    const char r2_107_zeta[] = "[-4140562374860211619063098135818149338786107 -43759983173909127950618820371560491745403 804423881624343280470169093105074305662 14435785514109985045742059041428530499 -12694304118831251292805753466971243 -1583684104370601660958596498596282 -9850000453768762162333255111918 18555806501448223479465201439 774753906602523228175267544 10061192331993210798407354 51736080696715567981487 -483514772866500635341 -8212929939596578078 -55238882790857992 -115556237092373 5357747611274 75239690974 52676849 -5232193 -25466 121 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: dense, rank 2, p = 107" <<endl;
    if( test_Fp(r2_107, r2_107_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // dense, rank 2, p = 211
    const char r2_211[] = "223 \n[ [0 3 0][3 1 0][0 0 0][2 0 0][0 2 2][0 2 1][0 1 0][4 0 0][1 0 2][0 1 1][1 1 2][1 1 0][1 3 0][2 2 0][2 0 1][1 0 1][1 0 0][2 0 2][0 2 0][0 1 2][1 2 1][3 0 0][0 0 3][0 3 1][0 4 0][1 0 3][0 0 1][1 1 1][0 1 3][0 0 2][3 0 1][1 2 0][2 1 0][2 1 1][0 0 4] ]\n[48 44 174 166 64 51 21 22 98 174 5 191 57 75 107 97 151 196 38 198 48 167 136 49 21 183 115 159 201 66 207 51 27 43 192] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n";
    const char r2_211_zeta[] = "[20625387851152557293083391544660009040933060233823 14516450658375183600271847494683189214194074045 -684439954168106030096097357853095916825550192 -483738513133044328153315055998727373463174 4450663604672357743176472217709027171959 -11740078091987226966965107406249082490 87994179164798562148297089463878453 529496229691210542228455906091331 -2646649344376391559509290896175 -10444176784983069831247425958 52917324575943572029498926 237297419623065345423762 -941802036543272594074 -4799235816058258225 19307623910855773 64522350049131 -173108140870 1319658473 -2884282 -82064 35 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: dense, rank 2, p = 211" <<endl;
    if( test_Fp(r2_211, r2_211_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    // dense, rank 2, p = 1009
    const char r2_1009[] = "1009 \n[ [0 3 0][3 1 0][0 0 0][2 0 0][0 2 2][0 2 1][0 1 0][4 0 0][1 0 2][0 1 1][1 1 2][1 1 0][1 3 0][2 2 0][2 0 1][1 0 1][1 0 0][2 0 2][0 2 0][0 1 2][1 2 1][3 0 0][0 0 3][0 3 1][0 4 0][1 0 3][0 0 1][1 1 1][0 1 3][0 0 2][3 0 1][1 2 0][2 1 0][2 1 1][0 0 4] ]\n[48 44 960 952 64 51 21 22 884 960 5 977 57 75 893 883 937 982 38 984 48 953 922 49 21 969 901 945 987 852 993 51 27 43 978] \n[ [-1 -1 -1][ 0  0  1][ 0  1  0][ 1  0  0] ] \n[4 0 0 0] \n";
    const char r2_1009_zeta[] = "[-1207020068576253752077404803423021173724331179038253595707180209 735061790287096337411258021829572625075102404429232280475180 598079304415313845938697421770600828292507950788834730589 758107537417407017989551836171239082478562633671283419 -298922374803585575685640183918239564901977018082819 93795358123965458388872713460307595815123791218 -289078976064337932862318909477851010815545055 -785348352159864027947534610952478114848971 -161458902868717003410496544813572449745 116979212631721420057500657873884554 519523589581200090948564920340475 -514889583331219118878657007275 -113876784378716241408664426 154385378152419311145505 737605545460107983259 266683417330332495 -84992091650002 266056125811 -662770731 -513581 -620 1]";
    user_time = get_cpu_time();
    get_timestamp(&wtime1);
    cout << "Testing: dense, rank 2, p = 1009" <<endl;
    if( test_Fp(r2_1009, r2_1009_zeta, 0))
        cout <<"PASS"<<endl;
    else
        val = false;
    get_timestamp(&wtime2);
    wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
    user_time = get_cpu_time() - user_time;
    printf("Time: CPU %.2f s, Wall: %.2f s\n\n\n", user_time, wall_time );

    return not val;
}

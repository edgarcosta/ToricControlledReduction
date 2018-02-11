
// Copyright 2018 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#define verbose 0

int main()
{

int examples_length = 10;

char examples[10][3][buffer_length] = {
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 127 = 2^7 - 1",
 "127 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[-4579937329576774398276408998492161 61902556751673186113476170975962 -67076381656801111876193207010 -2552775679852010819435174545 -1663497592083851742233076 183377687316330113317032 -103137056983312774644 -10489661307357926555 0 650360301776795 396460432884 -43704299688 24580596 2338705 3810 -218 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 251 = 2^8 - 5",
 "251 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[248187293863347975401129272743132504001 -760307736633167080719638571442113193 -3186063425810394017424326115284953 17008012321559617151873622816272 -1992970743093463458152522008 -202473123302323976824260204 537777219926491306306136 1937737229053152996123 -21842145425144859348 30757245584247123 135490085170136 -809701704204 -126506008 17136272 -50953 -193 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 509 = 2^9 - 3",
 "509 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[20299461509259581849945602588952505776043841 -20763225786351716992892511168601379609665 -61727052973473985351595040060860386429 92541133462429581829258299529287186 108729325883281414761446054764947 226453769003366661982893003794 -1350751118455258581693678521 -914614951405824170312363 6408598186512883363556 -3530227810629973523 -20123531898352361 13021855124834 24132617907 79278786 -204109 -265 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 1021 = 2^10 - 3",
 "1021 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[-1394478662688719756722453406066490116670622144321 1668118284270125154932412863044443930628463207 906652245675557955911605837341062389103012 -1246028956297756895922858792040801520411 -516565009202603295999555133674573831 321290532400798464021685456720461 -95248637087308583567733350899 1180880460795553919187361000 0 -1132803161805372121000 87650783332638979 -283624325243541 437440559271 1012210211 -706532 -1247 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 2039 = 2^11 - 9",
 "2039 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[89263849625258546797418573345596138821802616271519361 9489939205205284034514560339864427229408283540122 -25482339711149316122406329984390265274176239980 -2168982258902377354685017893523383163670820 1671604066974437812760820908621678329740 790004151318785537262682314476338236 -660974599230787619354976011183515 -173286926840823670781912198980 269611494298350741775649360 -41680349140948096421380 -38239822443328105915 10993247830420476 5594942310540 -1746158820 -4934380 442 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 4093 = 2^12 - 3",
 "4093 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[-6203944647560564735898307646931044205113365582336561108801 412173047542207119068442036303288644099877929343627537 391769458881050794887323023002092593594364852326690 2851613394786883511502185614445095218258236529 -5751866517135572625885805562122562074112205 -1305009694886627994664212390773486828261 430384842273967557838575550622596155 -335539425233563432295081644948260 0 20029036914315069353812740 -1533521922841506430155 277564084783511789 73025550860205 -2161091721 -17722690 -1113 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 8191 = 2^13 - 1",
 "8191 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[410573405487156078689423742750628108275098164385170350367703041 36533501134923327545091971546278782065047438644867924793170 313782939780655558117340212015947733362615875757833820 -740079056973170087542104423640755506480476602048434 -79217301832887157101381778791826030299794744374 -10662314152122536988201620439506974038843043 -753342581847243724903937166086656642629 171502739702380260603572311141916944 36849206910650013217496039540176 2556214007086431348448305424 -167357358457962643115589 -35304488095825436323 -3909530797450294 -544388390834 3440220 5970 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 16381 = 2^14 - 3",
 "16381 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[26881070824087011214685642408694287452112909430103080458284886710081 -746314736628849678965665273416127561691841761917828255678021450 -144562095463924505441271676248623605354160491303974625505099 3974398275082333554054894905214176095068989886516426806 422800263352070448953792846596373249763122334861912 -13342054500906058975901523256281551514920271590 -989100630378527138681606318382508620759366 38045297982471766718576636230675636458 1144488133348170707699163304252576 141781696730672971227330813978 -13736586882696960014963046 -690526338626247842390 81547746949094232 2856717416006 -387230459 -7450 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 32749 = 2^15 - 19",
 "32749 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[-1750526507940367569923766935249595258387043187585129894976283316831976001 -54360324823839709661165810363406107165963509638310635380681233128305 1554247975383896012166195609468896941254906862659622245075107565 72442360834832236890588189124010536058994717630422713229601 624053968488323901021209094696646884751507798897471321 -27532732539398260755623911423109131911434566474403 -1725375570994243425163453671262434254441905180 -27594044730690486556208018503990213748856 0 25728784980248617549475104316856 1500000734851890274632985180 22318297202238945601403 -471669510353323321 -51051929744601 -1021277565 33305 1]"
},
// f = 3*x + y + z + x^-2*y^2*z - 2 - x^-1*y - y^-1*z^-1 - x^2*y^-4*z^-1 - x*y^-3*z^-1 + x^3*y^-6*z^-2 + 3*x^-2*y^-1*z^-2
{
"K3_nonWPS, p = 65521 = 2^16 - 15",
 "65521 \n[ [-1  1  0][-2  2  1][ 1  0  0][ 2 -4 -1][ 0  0  1][ 3 -6 -2][ 0 -1 -1][ 0  0  0][-2 -1 -2][ 0  1  0][ 1 -3 -1] ]\n[-1 1 3 -1 1 1 -1 -2 3 1 -1] \n[ [ 1  1 -1][-1 -1  2][-1  0 -1][-1 -1 -1][ 0 -1  1] ] \n[1 1 1 1 1] \n",
"[115368773083406531811195955330938200073304864163668460305155754678164266776321 -1437930321294255838015140787811468678079654846033723543198775498748903267 -22258635705960388660114501315781503583407237448497314770573454822509 727742680413032942373120039032817714361172663771737180912479455 -5025886924532502206962452024390889735135439158344769688205 -63282684486300572132775675099638390274983061538140199 2476341570629406842202110613347681795106325845032 -2656479219208293968297677859133254970426381 -465967550900305196341915003906872547326 -618792994998273509384009046535341 134365718796676797821146066472 -799837553702006091198919 -14796826322806107405 499082882523455 -3555759149 -53507 1]"
},

};


    
    return not run_examples(examples, examples_length, verbose);
}

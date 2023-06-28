#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "include/Mjday.h"
#include "include/Position.h"
#include "include/MeanObliquity.h"
#include "include/Frac.h"
#include "include/Matrix.h"
#include "include/R_x.h"
#include "include/R_y.h"
#include "include/R_z.h"
#include "include/LTC.h"
#include "include/PrecMatrix.h"
#include "include/PoleMatrix.h"
#include "include/Gmst.h"
#include "include/NutAngles.h"
#include "include/NutMatrix.h"
#include "include/EqnEquinox.h"
#include "include/Gast.h"
#include "include/Timediff.h"
#include "include/GHAMatrix.h"
#include "include/Norm.h"
#include "include/Doubler.h"
#include "include/IERS.h"
#include "include/Geodetic.h"
#include "include/Anglesdr.h"
#include "include/MeasUpdate.h"
#include "include/TimeUpdate.h"
#include "include/Legendre.h"
#include "include/Globals.h"
#include "include/AccelHarmonic.h"
#include "include/Accel.h"
#include "include/G_AccelHarmonic.h"
#include "include/VarEqn.h"
#include "include/AzElPa.h"
#include "include/Sign_.h"
#include "include/DEInteg.h"


using namespace std;

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

double **eopdata = nullptr;
double **Cnm = nullptr;
double **Snm = nullptr;
AuxParam auxParam;

int Mjday_01(){
    _assert(fabs(Mjday(2023,4,27,18,33,34) - 60061.7733101854) < pow(10, -11));
    return 0;
}

int Mjday_02(){
    _assert(fabs(Mjday(2023,4,27,0,0,0) - 60061 == 0));
    return 0;
}

int Mjday_03(){
    _assert(fabs(Mjday(0,0,0,0,0,0) +678987)==0);
    return 0;
}

int Mjday_04(){
    _assert(fabs(Mjday(1,1,1,1,1,1) +678589.957627315)< pow(10, -9));
    return 0;
}

int Position_01(){
    Matrix v(3,1);
    v = Position(1,2,3);
    _assert(fabs(v(0,0)-(1.0e+06 *-1.438078943440528)) < 1e-8 &&
            fabs(v(1,0)-(1.0e+06 *-2.239675255177838)) < 1e-8 &&
            fabs(v(2,0)-(1.0e+06 *5.776811079007394)) < 1e-8);
    return 0;
}

int Position_02(){
    Matrix v(3,1);
    v = Position(10.3,-0.3,13);
    _assert(fabs(v(0,0)-(1.0e+06 *-3.905876347704072)) < 1e-8 &&
            fabs(v(1,0)-(1.0e+06 *-4.679092129311511))< 1e-8 &&
            fabs(v(2,0)-(1.0e+06 *-1.872801712574203)) < 1e-8);
    return 0;
}

int Position_03(){
    Matrix v(3,1);
    v = Position(1,1,1);
    _assert(fabs(v(0,0)-(1.0e+06 *1.86637669553512)) < 1e-8 &&
            fabs(v(1,0)-(1.0e+06 *2.90670948274229)) < 1e-8 &&
            fabs(v(2,0)-(1.0e+06 *5.34376928735908)) < 1e-8);
    return 0;
}

int MOblq_01(){
    _assert(fabs(MeanObliquity(4.974611324) - 0.409413039277443) < pow(10, -9));
    return 0;
}

int MOblq_02(){
    _assert(fabs(MeanObliquity(0) - 0.409413070181314) < pow(10, -9));
    return 0;
}

int MOblq_03(){
    _assert(fabs(MeanObliquity(1.123456789) -0.409413063202043) < pow(10, -9));
    return 0;
}

int Frac_01(){
    _assert(Frac(2.13212313221312)- 0.132123132213120 < 1e-8);
    return 0;
}

int Frac_02(){
    _assert(Frac(1) == 0);
    return 0;
}

int Frac_03(){
    _assert(Frac(15.3912930312) - 0.3912930312 < 1e-8);
    return 0;
}

int R_x_01(){
    Matrix m(3,3);
    m = R_x(1.9);
    _assert(fabs(m(0,0)-1) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2) - 0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) +0.323289566863503) < 1e-8 && fabs(m(1,2)-0.946300087687414) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) +0.946300087687414) < 1e-8 && fabs(m(2,2)+0.323289566863503) < 1e-8);
    return 0;
}

int R_x_02(){
    Matrix m(3,3);
    m= R_x(1.0);
    _assert(fabs(m(0,0)-1) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2) - 0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -0.54030230586814) < 1e-8 && fabs(m(1,2)-0.841470984807897) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) +0.841470984807897) < 1e-8 && fabs(m(2,2)-0.54030230586814) < 1e-8);
    return 0;
}

int R_x_03(){
    Matrix m(3,3);
    m = R_x(-1.69420);
    _assert(fabs(m(0,0)-1) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2) - 0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) +0.123090703484145) < 1e-8 && fabs(m(1,2)+0.992395424574186) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0.992395424574186) < 1e-8 && fabs(m(2,2)+0.123090703484145) < 1e-8);
    return 0;
}

int R_x_04(){
    Matrix m(3,3);
    m = R_x(-7.5421);
    _assert(fabs(m(0,0)-1) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2) - 0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -0.306850038549201) < 1e-8 && fabs(m(1,2)+0.951757875639784) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0.951757875639784) < 1e-8 && fabs(m(2,2)-0.306850038549201) < 1e-8);
    return 0;
}

int R_y_01(){
    Matrix m(3,3);
    m = R_y(-1.9873);
    _assert(fabs(m(0,0)+0.404565509944959) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2) - 0.914509020274254) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -1.0) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)+0.914509020274254) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)+0.404565509944959) < 1e-8);
    return 0;
}

int R_y_02(){
    Matrix m(3,3);
    m = R_y(0.83);
    _assert(fabs(m(0,0)-0.674875760071267) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2)+0.737931371109963) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -1.0) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)-0.737931371109963) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)-0.674875760071267) < 1e-8);
    return 0;
}

int R_y_03(){
    Matrix m(3,3);
    m=R_y(0.0);
    _assert(fabs(m(0,0)-1) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -1.0) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)-1) < 1e-8);
    return 0;
}

int R_y_04(){
    Matrix m(3,3);
    m=R_y(556.894);
    _assert(fabs(m(0,0)+0.673324384626954) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2)-0.739347193858700) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -1.0) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)+0.739347193858700) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)+0.673324384626954) < 1e-8);
    return 0;
}

int R_z_01(){
    Matrix m(3,3);
    m = R_z(8.36);
    _assert(fabs(m(0,0)+0.484698437250152) < 1e-8 && fabs(m(0,1) - 0.874681327642965) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)+0.874681327642965) < 1e-8 && fabs(m(1,1) +0.484698437250152) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)-1) < 1e-8);
    return 0;
}

int R_z_02(){
    Matrix m(3,3);
    m = R_z(12.92);
    _assert(fabs(m(0,0)-0.938122020292311) < 1e-8 && fabs(m(0,1)-0.346304887408008) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)+0.346304887408008) < 1e-8 && fabs(m(1,1)-0.938122020292311) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)-1) < 1e-8);
    return 0;
}

int R_z_03(){
    Matrix m(3,3);
    m = R_z(0.0);
    _assert(fabs(m(0,0)-1) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -1.0) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)-1) < 1e-8);
    return 0;
}

int R_z_04(){
    Matrix m(3,3);
    m = R_z(-1.0);
    _assert(fabs(m(0,0)-0.540302305868140) < 1e-8 && fabs(m(0,1)+0.841470984807897) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)-0.841470984807897) < 1e-8 && fabs(m(1,1) -0.540302305868140) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)-1) < 1e-8);
    return 0;
}

int LTC_01(){
    Matrix m(3,3);
    m = LTC(0.56,34.34);
    _assert(fabs(m(0,0)+0.531186197920883) < 1e-8 && fabs(m(0,1)-0.847255111013416) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)+0.182844377938121) < 1e-8 && fabs(m(1,1) +0.114634197735305) < 1e-8 && fabs(m(1,2)+0.976435832078076) < 1e-8 &&
            fabs(m(2,0)+0.827290249304788) < 1e-8 && fabs(m(2,1) +0.518669237155267) < 1e-8 && fabs(m(2,2)-0.215807937374869) < 1e-8);
    return 0;
}

int LTC_02(){
    Matrix m(3,3);
    m = LTC(-12.69,420.69);
    _assert(fabs(m(0,0)-0.123314696335228) < 1e-8 && fabs(m(0,1)-0.992367616192583) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)-0.277502299096376) < 1e-8 && fabs(m(1,1) +0.034483301537680) < 1e-8 && fabs(m(1,2)-0.960105919110640) < 1e-8 &&
            fabs(m(2,0)-0.952778022240215) < 1e-8 && fabs(m(2,1) +0.118395169864784) < 1e-8 && fabs(m(2,2)+0.279636592899987) < 1e-8);
    return 0;
}

int LTC_03(){
    Matrix m(3,3);
    m = LTC(0.0,0.0);
    _assert(fabs(m(0,0)-0) < 1e-8 && fabs(m(0,1)-1) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1)-0) < 1e-8 && fabs(m(1,2)-1) < 1e-8 &&
            fabs(m(2,0)-1) < 1e-8 && fabs(m(2,1)-0) < 1e-8 && fabs(m(2,2)-0) < 1e-8);
    return 0;
}

int PrecMatrix_01(){
    Matrix m(3,3);
    m = PrecMatrix(51544.5, 67391.6);
    _assert(fabs(m(0,0)-0.999944037343594) < 1e-8 && fabs(m(0,1)+0.009703172634268) < 1e-8 && fabs(m(0,2)+0.004215521536129) < 1e-8 &&
            fabs(m(1,0)-0.009703172627839) < 1e-8 && fabs(m(1,1)-0.999952922903166) < 1e-8 && fabs(m(1,2)+0.000020454064126) < 1e-8 &&
            fabs(m(2,0)-0.004215521550929) < 1e-8 && fabs(m(2,1)+0.000020451013719) < 1e-8 && fabs(m(2,2)-0.999991114440428) < 1e-8);
    return 0;
}

int PrecMatrix_02(){
    Matrix m(3,3);
    m = PrecMatrix(15663.8, 54355.4);
    _assert(fabs(m(0,0)-0.999666608088714) < 1e-8 && fabs(m(0,1)+0.023679054486938) < 1e-8 && fabs(m(0,2)+0.010294418439651) < 1e-8 &&
            fabs(m(1,0)-0.023679054258412) < 1e-8 && fabs(m(1,1)-0.999719604451204) < 1e-8 && fabs(m(1,2)+0.000121923565862) < 1e-8 &&
            fabs(m(2,0)-0.010294418965303 ) < 1e-8 && fabs(m(2,1)+0.000121879175260) < 1e-8 && fabs(m(2,2)-0.999947003637510) < 1e-8);
    return 0;
}

int PrecMatrix_03(){
    Matrix m(3,3);
    m = PrecMatrix(78113.5, 46712.9);
    _assert(fabs(m(0,0)-0.999780267803598) < 1e-8 && fabs(m(0,1)-0.019226462983197) < 1e-8 && fabs(m(0,2)-0.008352199214632) < 1e-8 &&
            fabs(m(1,0)+0.019226463082308) < 1e-8 && fabs(m(1,1)-0.999815151251014) < 1e-8 && fabs(m(1,2)+0.000080288580804) < 1e-8 &&
            fabs(m(2,0)+0.008352198986483) < 1e-8 && fabs(m(2,1)+0.000080312311038) < 1e-8 && fabs(m(2,2)-0.999965116552584) < 1e-8);
    return 0;
}

int PoleMatrix_01(){
    Matrix m(3,3);
    m = PoleMatrix(0,0);
    _assert(fabs(m(0,0)-1) < 1e-8 && fabs(m(0,1) - 0) < 1e-8 && fabs(m(0,2)-0) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -1.0) < 1e-8 && fabs(m(1,2)-0) < 1e-8 &&
            fabs(m(2,0)-0) < 1e-8 && fabs(m(2,1) -0) < 1e-8 && fabs(m(2,2)-1) < 1e-8);
    return 0;
}

int PoleMatrix_02(){
    Matrix m(3,3);
    m=PoleMatrix(-1.93, 100.38);
    _assert(fabs(m(0,0)+0.351528841940960) < 1e-8 && fabs(m(0,1) -0.140793672822927) < 1e-8 && fabs(m(0,2)+0.925529370131860) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) -0.988626422578826) < 1e-8 && fabs(m(1,2)-0.150392142677043) < 1e-8 &&
            fabs(m(2,0)-0.936177052316306) < 1e-8 && fabs(m(2,1) -0.052867175752281) < 1e-8 && fabs(m(2,2)+0.347530701441369) < 1e-8);
    return 0;
}

int PoleMatrix_03(){
    Matrix m(3,3);
    m=PoleMatrix(69.431, -312.3123);
    _assert(fabs(m(0,0)-0.950498013319686) < 1e-8 && fabs(m(0,1) -0.298956146610814) < 1e-8 && fabs(m(0,2)+0.084727498953669) < 1e-8 &&
            fabs(m(1,0)-0) < 1e-8 && fabs(m(1,1) +0.272671857178060) < 1e-8 && fabs(m(1,2)+0.962107092949152) < 1e-8 &&
            fabs(m(2,0)+0.310730633628760) < 1e-8 && fabs(m(2,1) -0.914480880448947) < 1e-8 && fabs(m(2,2)+0.259174058535935) < 1e-8);
    return 0;
}

int Gmst_01(){
    _assert(fabs(gmst(0)-0.973208148169487) < 1e-9);
    return 0;
}

int Gmst_02(){
    _assert(fabs(gmst(6.008176636574091e+04)-2.636955396305363) < 1e-9);
    return 0;
}

int Gmst_03(){
    _assert(fabs(gmst(54133.32)-4.316622222700707) < 1e-9);
    return 0;
}

int NutAngles_01(){
    double dpsi;
    double desp;
    NutAngles(183.92, dpsi,desp);
    _assert(fabs(dpsi-42.2232973928201e-006)<1e-9 && fabs(desp - 34.5749583506206e-006) < 1e-9);
    return 0;
}

int NutAngles_02(){
    double dpsi;
    double desp;
    NutAngles(51447.87, dpsi,desp);
    _assert(fabs(dpsi+6.705345712767994e-05)<1e-9 && fabs(desp+2.559670567237168e-05)<1e-9);
    return 0;
}

int NutAngles_03(){
    double dpsi;
    double desp;
    NutAngles(-73621.93, dpsi, desp);
    _assert(fabs(dpsi-7.873246165551288e-05)<1e-9 && fabs(desp-4.127632047079497e-06)<1e-9);
    return 0;
}

int NutMatrix_01(){
    Matrix m(3,3);
    m = NutMatrix(5144.87);
    _assert(fabs(m(0,0)-0.999999997347533) < 1e-8 && fabs(m(0,1) -0.000066816427691) < 1e-8 && fabs(m(0,2)-0.000028991353905) < 1e-8 &&
            fabs(m(1,0)+0.000066815904891) < 1e-8 && fabs(m(1,1)-0.999999997605213) < 1e-8 && fabs(m(1,2)+0.000018033543150) < 1e-8 &&
            fabs(m(2,0)+0.000028992558772) < 1e-8 && fabs(m(2,1) -0.000018031606019) < 1e-8 && fabs(m(2,2)-0.999999999417146) < 1e-8);
    return 0;
}

int NutMatrix_02(){
    Matrix m(3,3);
    m = NutMatrix(0.0);
    _assert(fabs(m(0,0)-0.999999999612860) < 1e-8 && fabs(m(0,1)+0.000025526217598) < 1e-8 && fabs(m(0,2)+0.000011076682963) < 1e-8 &&
            fabs(m(1,0)-0.000025525777883) < 1e-8 && fabs(m(1,1)-0.999999998886346) < 1e-8 && fabs(m(1,2)+0.000039695628802) < 1e-8 &&
            fabs(m(2,0)-0.000011077696230) < 1e-8 && fabs(m(2,1)-0.000039695346046) < 1e-8 && fabs(m(2,2)-0.999999999150782) < 1e-8);
    return 0;
}

int NutMatrix_03(){
    Matrix m(3,3);
    m = NutMatrix(-8311.932);
    _assert(fabs(m(0,0)-0.999999998221923) < 1e-8 && fabs(m(0,1)-0.000054703835413) < 1e-8 && fabs(m(0,2)-0.000023741187621) < 1e-8 &&
            fabs(m(1,0)+0.000054703209493) < 1e-8 && fabs(m(1,1)-0.999999998156245) < 1e-8 && fabs(m(1,2)+0.000026364165549) < 1e-8 &&
            fabs(m(2,0)+0.000023742629798) < 1e-8 && fabs(m(2,1)-0.000026362866783) < 1e-8 && fabs(m(2,2)-0.999999999370643) < 1e-8);
    return 0;
}

int EqnEquinox_01(){
    _assert(fabs(EqnEquinox(0.0)-2.552621760089941e-05) < 1e-9);
    return 0;
}

int EqnEquinox_02(){
    _assert(fabs(EqnEquinox(18231.93)+8.046522659905191e-05) < 1e-9);
    return 0;
}

int EqnEquinox_03(){
    _assert(fabs(EqnEquinox(-51147.69)+2.884121114458045e-05) < 1e-9);
    return 0;
}

int Gast_01(){
    _assert(fabs(gast(28.9382132)-1.08283939664680) < 1e-5); //Corresponds to 18/05/2023 a las 18:06:23
    return 0;
}

int Gast_02(){
    _assert(fabs(gast(-12.93213)-1.17720356850459)<1e-5);
    return 0;
}

int Gast_03(){
    _assert(fabs(gast(0.0)-973.233674387088e-003)<1e-5);
    return 0;
}

int Timediff_01(){
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    timediff(15677.32, -13233.78, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    _assert(fabs(UT1_TAI-2.89111e+04)<1e-9 && fabs(UTC_GPS-1.325278e+04)<1e-9 && fabs(UT1_GPS-2.89301e+04)<1e-9 &&
                    fabs(TT_UTC+1.3201596e+04)<1e-9 && fabs(GPS_UTC+1.325278e+04)<1e-9);
    return 0;
}

int Timediff_02(){
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    timediff(-1983.63, 2034.67, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    _assert(fabs(UT1_TAI+4.0183e+03)<1e-9 && fabs(UTC_GPS+2.01567e+03)<1e-9 && fabs(UT1_GPS+3.9993e+03)<1e-9 &&
            fabs(TT_UTC-2.066854e+03)<1e-9 && fabs(GPS_UTC-2.01567e+03)<1e-9);
    return 0;
}

int Timediff_03(){
    double UT1_TAI;
    double UTC_GPS;
    double UT1_GPS;
    double TT_UTC;
    double GPS_UTC;
    timediff(-2019.73112, -9381.98321, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    _assert(fabs(UT1_TAI-7.36225209e+03)<1e-9 && fabs(UTC_GPS-9.40098321e+03)<1e-9 && fabs(UT1_GPS-7.38125209e+03)<1e-9 &&
            fabs(TT_UTC+9.34979921e+03)<1e-9 && fabs(GPS_UTC+9.40098321e+03)<1e-9);
    return 0;
}

int GHAMatrix_01(){
    Matrix m(3,3);
    m = GHAMatrix(0);
    _assert(fabs(m(0,0)-0.562629168438446) < 1e-5 && fabs(m(0,1)-0.826709392000758) < 1e-5 && fabs(m(0,2)-0.0) < 1e-5 &&
            fabs(m(1,0)+0.826709392000758) < 1e-5 && fabs(m(1,1)-0.562629168438446) < 1e-5 && fabs(m(1,2)-0.0) < 1e-5 &&
            fabs(m(2,0)-0.0) < 1e-5 && fabs(m(2,1)-0.0) < 1e-5 && fabs(m(2,2)-1) < 1e-5);
    return 0;
}

int Norm_01() {
    Matrix v1 = Matrix(3, 1);
    v1(0,0) = 1;
    v1(1,0) = 1;
    v1(2,0) = 1;
    _assert(fabs(norm(v1) - sqrt(3)) < 1e-9);
    return 0;
}

int Norm_02(){
    double data[] = {-2,5,-3};
    Matrix v = Matrix(3,1,data,3);
    _assert(fabs(norm(v)-sqrt(38))<1e-9);
    return 0;
}

int Norm_03(){
    double data[] = {-1.5,2.7,-3.8};
    Matrix v = Matrix(3,1,data,3);
    _assert(fabs(norm(v)-sqrt(23.98))<1e-9);
    return 0;
}

int Dot_01(){
    double w[] = {4, -1, 2};
    double u[] = {2, -2, -1};
    Matrix m1 = Matrix(3,1,w,3);
    Matrix m2 = Matrix(3,1,u,3);
    _assert(m1.dot(m2) == 8);
    return 0;
}

int Dot_02(){
    double w[] = {1, -6, 9};
    double u[] = {0.232, -5.3, -5.6};
    Matrix m1 = Matrix(3,1,w,3);
    Matrix m2 = Matrix(3,1,u,3);
    _assert(fabs(m1.dot(m2)+18.368)<1e-9);
    return 0;
}

int Dot_03(){
    double w[] = {0, 19.321, -9.83};
    double u[] = {-1233.232, 50.3312, 9.312};
    Matrix m1 = Matrix(3,1,w,3);
    Matrix m2 = Matrix(3,1,u,3);
    _assert(fabs(m1.dot(m2)-880.9121552)<1e-9);
    return 0;
}

int Transpose_01(){
    double data[] = {1.5,3.4,-5.6,9.8,-1,2,5.4,65.13,3};
    Matrix m1 = Matrix(3,3,data,9);
    Matrix t = m1.transpose();
    _assert(t(0,0)==1.5 && t(0,1)==9.8 && t(0,2)==5.4 &&
            t(1,0)==3.4 && t(1,1)==-1 && t(1,2)==65.13 &&
            t(2,0)==-5.6 && t(2,1)==2 && t(2,2)==3);
    return 0;
}

int Transpose_02(){
    double data[] = {1.5,3.4,-5.6};
    Matrix m1 = Matrix(3,1,data,9);
    Matrix t = m1.transpose();
    _assert(t(0,0)==1.5 && t(0,1)==3.4 && t(0,2)==-5.6);
    return 0;
}

int Cross_01() {
    Matrix m1 = Matrix(3, 1);
    m1(0,0) = 4;
    m1(1,0) = -2;
    m1(2,0) = 1;
    Matrix m2 = Matrix(3, 1);
    m2(0,0) = 1;
    m2(1,0) = -1;
    m2(2,0) = 3;
    Matrix res = m1.cross(m2);
    _assert(res(0, 0) == -5 && res(1, 0) == -11 && res(2, 0) == -2);
    return 0;
}


int Cross_02(){
    double data1[] = {0.5, -1.2, 2.8};
    double data2[] = {-0.8, 1.6, -3.5};
    Matrix m1 = Matrix(3,1, data1, 3);
    Matrix m2 = Matrix(3,1, data2, 3);
    Matrix res = m1.cross(m2);
    _assert(fabs(res(0, 0)+0.28)<1e-9 && fabs(res(1, 0)+0.49)<1e-9 && fabs(res(2, 0)+0.16)<1e-9);
    return 0;
}

int Cross_03(){
    double data1[] = {67, -18.3, 2.6};
    double data2[] = {-10.2, 13.6, -4.5};
    Matrix m1 = Matrix(3,1, data1, 3);
    Matrix m2 = Matrix(3,1, data2, 3);
    Matrix res = m1.cross(m2);
    _assert(fabs(res(0, 0)-46.9900)<1e-9 && fabs(res(1, 0)-274.9800)<1e-9 && fabs(res(2, 0)-724.5400)<1e-9);
    return 0;
}

int Doubler_01(){
    double data1[] = {-0.0514407301715203, 0.838593164440367, 0.54232403213698};
    double data2[] = {0.185350425424354, 0.924321659182723, 0.333578611665541};
    double data3[] = {0.48999206372453, 0.865773547227108, -0.10170517395279};
    double data4[] = {5854667.7577933, 962016.736146505, 2333503.53479825};
    double data5[] = {5847642.87233096, 1003838.42368066, 2333501.82312028};
    double data6[] = {5839555.2146941, 1049868.17436044, 2333499.77773523};
    double cc1=3542174.64126966;
    double cc2=5580277.97983915;
    double magrsite1=6375566.60240377;
    double magrsite2=6375566.60240377;
    double magr1in=7232409.73786304;
    double magr2in=7230730.23983205;
    Matrix los1(3,1, data1, 3);
    Matrix los2(3,1, data2, 3);
    Matrix los3(3,1, data3, 3);
    Matrix rsite1(3,1, data4, 3);
    Matrix rsite2(3,1, data5, 3);
    Matrix rsite3(3,1, data6, 3);
    double t1=-97.9999914765358;
    double t3=108.000017702579;
    char direct='y';

    Matrix r2(3,1);
    Matrix r3(3,1);
    double f1;
    double f2;
    double q1;
    double magr1;
    double magr2;
    double a;
    double deltae32;

    doubler(cc1,cc2,magrsite1,magrsite2,magr1in,magr2in,los1,los2,los3,rsite1,rsite2,
            rsite3,t1,t3,direct,r2,r3,f1,f2,q1,magr1,magr2,a,deltae32);

    _assert(fabs(r2(0,0) -6147304.28873136) < pow(10, -6) &&
            fabs(r2(1,0) -2498216.09757119) < pow(10, -6) &&
            fabs(r2(2,0) -2872808.05359544) < pow(10, -6) &&
            fabs(r3(0,0) -6515290.32938719) < pow(10, -6) &&
            fabs(r3(1,0) -2243833.6087938) < pow(10, -6) &&
            fabs(r3(2,0) -2193240.8583081) < pow(10, -6) &&
            fabs(f1 +2.05106789885611817680910462514e-009) < pow(10, -6) &&
            fabs(f2 -4.62715003743396913193919317564e-008) < pow(10, -6) &&
            fabs(q1 -4.6316936712375920737576871257e-008) < pow(10, -6) &&
            fabs(magr1 -7232409.73786304) < pow(10, -6) &&
            fabs(magr2 -7230730.23983205) < pow(10, -6) &&
            fabs(a -7458443.29599821) < pow(10, -4) &&
            fabs(deltae32 -0.10918876174611) < pow(10, -6));
    return 0;
}

int Doubler_02() {
    double data1[] = {-0.123456, 0.987654, 0.246813};
    double data2[] = {0.135791, 0.579135, 0.975312};
    double data3[] = {0.246813, 0.975312, -0.579135};
    double data4[] = {1234567.89, 987654.32, 54321.67};
    double data5[] = {987654.32, 54321.67, 1234567.89};
    double data6[] = {54321.67, 1234567.89, 987654.32};
    double cc1 = 1234567.89;
    double cc2 = 987654.32;
    double magrsite1 = 6378.16;
    double magrsite2 = 6378.16;
    double magr1in = 12345.67;
    double magr2in = 76543.21;
    Matrix los1(3, 1, data1, 3);
    Matrix los2(3, 1, data2, 3);
    Matrix los3(3, 1, data3, 3);
    Matrix rsite1(3, 1, data4, 3);
    Matrix rsite2(3, 1, data5, 3);
    Matrix rsite3(3, 1, data6, 3);
    double t1 = -123.456;
    double t3 = 789.012;
    char direct = 'x';

    Matrix r2(3, 1);
    Matrix r3(3, 1);
    double f1;
    double f2;
    double q1;
    double magr1;
    double magr2;
    double a;
    double deltae32;

    doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3, rsite1, rsite2,
            rsite3, t1, t3, direct, r2, r3, f1, f2, q1, magr1, magr2, a, deltae32);

    _assert(fabs(r2(0, 0) - 9.884495373199758e+05) < pow(10, -6) &&
            fabs(r2(1, 0) - 5.771319213772797e+04) < pow(10, -6) &&
            fabs(r2(2, 0) - 1.240279498241933e+06) < pow(10, -6) &&
            fabs(r3(0, 0) - -1.009390922232133e+06) < pow(10, -6) &&
            fabs(r3(1, 0) - -2.968823567318317e+06) < pow(10, -6) &&
            fabs(r3(2, 0) - 3.483605473717008e+06) < pow(10, -6) &&
            fabs(f1 + 67.467620416185412) < pow(10, -6) &&
            fabs(f2 - 5.988464397414708e+02) < pow(10, -6) &&
            fabs(q1 - 6.026349958272068e+02) < pow(10, -6) &&
            fabs(magr1 - 1.581999108305839e+06) < pow(10, -6) &&
            fabs(magr2 - 1.587027578252521e+06) < pow(10, -6) &&
            fabs(a - -1.805141938468124e+06) < pow(10, -6) &&
            fabs(deltae32 - 1.030790802175237) < pow(10, -6));

    return 0;
}

int IERS_01()
{

    double UT1_UTC;
    double TAI_UTC;
    double x_pole;
    double y_pole;

    IERS(eopdata, 49846, UT1_UTC, TAI_UTC, x_pole, y_pole);

    _assert(fabs(UT1_UTC-0.0538074)< 1e-7&&
            (fabs(TAI_UTC-29))< 1e-7&&
            (fabs(x_pole-0.000000793552730))< 1e-7&&
            (fabs(y_pole-0.000002583974502))< 1e-7);
    return 0;
}

int Geodetic_01(){
    double data[] = {1.4, -12.3, 9};
    Matrix r = Matrix(3,1,data,3);
    double lat, lon, h;
    Geodetic(r, lat, lon, h);
    _assert(fabs(lat+1.45746293002949)<1e-8 &&
            fabs(lon-1.57050742764477)<1e-8 &&
            fabs(h+6.35674331245698e06)<1e-8);
    return 0;
}

int Geodetic_02(){
    double data[] = {-10.69,-42.3,-73.4};
    Matrix r = Matrix(3,1,data,3);
    double lat, lon, h;
    Geodetic(r, lat, lon, h);
    _assert(fabs(lat+1.81833210258717)<1e-8 &&
            fabs(lon+1.56977966182565)<1e-8 &&
            fabs(h+6.35667889206670e+006)<1e-8);
    return 0;
}

int Geodetic_03(){
    double data[] = {65.9,-3.9,10};
    Matrix r = Matrix(3,1,data,3);
    double lat, lon, h;
    Geodetic(r, lat, lon, h);
    _assert(fabs(lat+59.1116312756943e-003)<1e-8 &&
            fabs(lon-1.56925575952758)<1e-8 &&
            fabs(h+6.35674226339469e+006)<1e-8);
    return 0;
}

int Anglesdr_01(){
    double obs12 =1.055908489493301294359639541653;
    double obs92 =1.363102145807571385915935024968;
    double obs182=1.976156026887587513485300405591;
    double obs13=0.282624656433945797839868419032;
    double obs93 =0.453434794338874846975073751310;
    double obs183=0.586427138011590742827650046820;
    double Mjd1 =49746.110150462947785854339599609375;
    double Mjd2 =49746.111284722108393907546997070312;
    double Mjd3 =49746.112534722313284873962402343750;

    double dataRs[] = {-5512568.445011530071496963500976562500, -2196994.687777972314506769180297851562, 2330805.221940449904650449752807617188};
    Matrix Rs = Matrix(3,1, dataRs, 3);

    Matrix r2 = Matrix(3,1);
    Matrix v2 = Matrix(3,1);

    anglesdr(obs12, obs92, obs182, obs13, obs93, obs183, Mjd1, Mjd2, Mjd3, Rs, Rs, Rs, r2, v2);

    _assert(fabs(r2(0,0)-6147304.288731360808014869689941406250)<pow(10, -4) &&
            fabs(r2(1,0)-2498216.097571192774921655654907226562)<pow(10, -4) &&
            fabs(r2(2,0)-2872808.053595441859215497970581054688)<pow(10, -4) &&
            fabs(v2(0,0)-3764.628994742533450335031375288963)<pow(10, -4) &&
            fabs(v2(1,0)+2217.844940728065466828411445021629)<pow(10, -4) &&
            fabs(v2(2,0)+6141.471007388881844235584139823914)<pow(10, -4));
    return 0;
}

int TimeUpdate_01(){
    double data1[] = {1,2,3,4,5,6,7,8,9};
    double data2[] = {0.1, 0.2, 0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    double data3[] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
    Matrix P = Matrix(3,3,data1,9);
    Matrix Phi = Matrix(3,3,data2,9);
    Matrix Qdt = Matrix(3,3,data3,9);
    P = TimeUpdate(P, Phi, Qdt);
    _assert(fabs(P(0,0)-2.29) < 1e-9 && fabs(P(0,1) -5.539999999999999) < 1e-9 && fabs(P(0,2)-8.789999999999999) < 1e-9 &&
            fabs(P(1,0)-5.2) < 1e-9 && fabs(P(1,1) -12.5) < 1e-9 && fabs(P(1,2)-19.800000000000001) < 1e-9 &&
            fabs(P(2,0)-8.11) < 1e-9 && fabs(P(2,1) -19.460000000000001) < 1e-9 && fabs(P(2,2)-30.809999999999999) < 1e-9);
    return 0;
}

int TimeUpdate_02(){
    double data1[] = {2, 4, 6,8, 10, 12, 14, 16, 18};
    double data2[] = {0.5, 0.3, 0.1,0.7, 0.9, 0.2, 0.4, 0.6, 0.8};
    Matrix P = Matrix(3,3,data1,9);
    Matrix Phi = Matrix(3,3,data2,9);
    Matrix Qdt = Matrix(3,3);
    P = TimeUpdate(P, Phi, Qdt);
    _assert(fabs(P(0,0)-5.22) < 1e-9 && fabs(P(0,1) -10.979999999999999) < 1e-9 && fabs(P(0,2)-12.6) < 1e-9 &&
            fabs(P(1,0)-12.059999999999999) < 1e-9 && fabs(P(1,1) -25.200000000000003) < 1e-9 && fabs(P(1,2)-28.440000000000005) < 1e-9 &&
            fabs(P(2,0)-16.920000000000002) < 1e-9 && fabs(P(2,1) -34.920000000000002) < 1e-9 && fabs(P(2,2)-38.160000000000004) < 1e-9);
    return 0;
}

int TimeUpdate_03(){
    double data1[] = {1,0,0,0,2,0,0,0,3};
    double data2[] = {0.8, 0.2, 0.1,0.3, 0.9, 0.4, 0.5, 0.7, 0.6};
    double data3[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6,  0.7, 0.8, 0.9};
    Matrix P = Matrix(3,3,data1,9);
    Matrix Phi = Matrix(3,3,data2,9);
    Matrix Qdt = Matrix(3,3,data3,9);
    P = TimeUpdate(P, Phi, Qdt);
    _assert(fabs(P(0,0)-0.85) < 1e-9 && fabs(P(0,1) -0.92) < 1e-9 && fabs(P(0,2)-1.16) < 1e-9 &&
            fabs(P(1,0)-1.12) < 1e-9 && fabs(P(1,1) -2.69) < 1e-9 && fabs(P(1,2)-2.73) < 1e-9 &&
            fabs(P(2,0)-1.56) < 1e-9 && fabs(P(2,1) -2.93) < 1e-9 && fabs(P(2,2)-3.21) < 1e-9);
    return 0;
}

int MeasUpdate_01(){
    double data1[] = {1,2,3};
    double data2[] = {4,5,6};
    double data3[] = {7,8,9};
    double data4[] = {0.1,0.2,0.3};
    Matrix x = Matrix(3,1,data1,3);
    Matrix z = Matrix(3,1,data2,3);
    Matrix g = Matrix(3,1,data3,3);
    Matrix s = Matrix(3,1,data4,3);
    Matrix G = Matrix::identity(3);
    Matrix P = Matrix::identity(3);
    Matrix K = MeasUpdate(x,z,g,s,G,P,3);
    _assert(fabs(K(0,0)-0.990099009900990)<1e-9 && fabs(K(1,1)-0.961538461538461)<1e-9 && fabs(K(2,2)-0.917431192660550)<1e-9
            && fabs(x(0,0)+1.970297029702970)<1e-9 && fabs(x(1,0)+0.884615384615384)<1e-9 && fabs(x(2,0)-0.247706422018349)<1e-9
            && fabs(P(0,0)-0.009900990099010)<1e-9 && fabs(P(1,1)-0.038461538461539)<1e-9 && fabs(P(2,2)-0.082568807339450)<1e-9);
    return 0;
}

int MeasUpdate_02(){
    double data1[] = {-10,28,23};
    double data2[] = {4.3,33,-12.2};
    double data3[] = {-13,44,-9};
    double data4[] = {0.5,0.6,0.7};
    Matrix x = Matrix(3,1,data1,3);
    Matrix z = Matrix(3,1,data2,3);
    Matrix g = Matrix(3,1,data3,3);
    Matrix s = Matrix(3,1,data4,3);
    Matrix G = Matrix::identity(3);
    Matrix P = Matrix::identity(3);
    Matrix K = MeasUpdate(x,z,g,s,G,P,3);
    _assert(fabs(K(0,0)-0.8)<1e-9 && fabs(K(1,1)-0.735294117647059)<1e-9 && fabs(K(2,2)-0.671140939597315)<1e-9
            && fabs(x(0,0)-3.840000000000002)<1e-9 && fabs(x(1,0)-19.911764705882351)<1e-9 && fabs(x(2,0)-20.852348993288590)<1e-9
            && fabs(P(0,0)-0.2)<1e-9 && fabs(P(1,1)-0.264705882352941)<1e-9 && fabs(P(2,2)-0.328859060402685)<1e-9);
    return 0;
}

int MeasUpdate_03(){
    double data1[] = {83,-2.8,13.65};
    double data2[] = {-10.3,3.33,-0.45};
    double data3[] = {-1.36,4.47,-0.91};
    double data4[] = {0.1,-0.3,0.6};
    Matrix x = Matrix(3,1,data1,3);
    Matrix z = Matrix(3,1,data2,3);
    Matrix g = Matrix(3,1,data3,3);
    Matrix s = Matrix(3,1,data4,3);
    Matrix G = Matrix::identity(3);
    Matrix P = Matrix::identity(3);
    Matrix K = MeasUpdate(x,z,g,s,G,P,3);
    _assert(fabs(K(0,0)-0.990099009900990)<1e-9 && fabs(K(1,1)-0.917431192660550)<1e-9 && fabs(K(2,2)-0.735294117647059)<1e-9
            && fabs(x(0,0)-74.148514851485146)<1e-9 && fabs(x(1,0)+3.845871559633027)<1e-9 && fabs(x(2,0)-13.988235294117647)<1e-9
            && fabs(P(0,0)-0.009900990099010)<1e-9 && fabs(P(1,1)-0.082568807339450)<1e-9 && fabs(P(2,2)-0.264705882352941)<1e-9);
    return 0;
}

int MeasUpdate_04(){
    Matrix x(6,1);
    x(0,0)=5.74789643764892e+006;
    x(1,0)=2.70256057450949e+006;
    x(2,0)=3.45908529959328e+006;
    x(3,0)=4.37955774313217e+003;
    x(4,0)=-1.94898340947690e+003;
    x(5,0)=-5.81331374009477e+003;

    Matrix z(1,1);
    z(0,0)=1.05590848949330e+000;

    Matrix g(1,1);
    g(0,0)=1.05595496938417e+000;

    Matrix s(1,1);
    g(0,0)=390.953752446730e-006;

    Matrix G(1,6);
    G(0,0)=118.226811703679e-009;
    G(0,1)=268.617065403109e-009;
    G(0,2)=-404.487218792015e-009;
    G(0,3)=0.00000000000000e+000;
    G(0,4)=0.00000000000000e+000;
    G(0,5)=0.00000000000000e+000;

    Matrix P(6,6);
    P(0,0)=101.488000600326e+006;
    P(1,0)=128.617803934412e+003;
    P(2,0)=169.139099218971e+003;
    P(3,0)=40.3080399813276e+003;
    P(4,0)=3.49644199102555e+003;
    P(5,0)=4.57197846555939e+003;
    P(0,1)=128.617803934412e+003;
    P(1,1)=101.286995428481e+006;
    P(2,1)=82.5095941203040e+003;
    P(3,1)=3.49657950927422e+003;
    P(4,1)=34.7617247204656e+003;
    P(5,1)=2.21017635089282e+003;
    P(0,2)=169.139099218971e+003;
    P(1,2)=82.5095941203040e+003;
    P(2,2)=101.332421164292e+006;
    P(3,2)=4.57226963026297e+003;
    P(4,2)=2.21023022307775e+003;
    P(5,2)=35.9528041664988e+003;
    P(0,3)=40.3080399813276e+003;
    P(1,3)=3.49657950927422e+003;
    P(2,3)=4.57226963026297e+003;
    P(3,3)=1.00166889789894e+003;
    P(4,3)=1.41783559164539e+000;
    P(5,3)=1.84417805416141e+000;
    P(0,4)=3.49644199102555e+003;
    P(1,4)=34.7617247204656e+003;
    P(2,4)=2.21023022307775e+003;
    P(3,4)=1.41783559164539e+000;
    P(4,4)=999.390221437519e+000;
    P(5,4)=884.050439731921e-003;
    P(0,5)=4.57197846555939e+003;
    P(1,5)=2.21017635089282e+003;
    P(2,5)=35.9528041664988e+003;
    P(3,5)=1.84417805416141e+000;
    P(4,5)=884.050439731921e-003;
    P(5,5)=999.857769260102e+000;

    int n=6;

    Matrix K(6,1);

    K = MeasUpdate(x,z,g,s,G,P,n);

    _assert(fabs(K(0,0)-(470.444610235622e+003)) < pow(10, -9) &&
            fabs(K(1,0)-(1.06906109631921e+006)) < pow(10, -9) &&
            fabs(K(2,0)+(1.60994712085869e+006)) < pow(10, -9) &&
            fabs(K(3,0)-(151.587833026057)) < pow(10, -9) &&
            fabs(K(4,0)-(348.248949621027)) < pow(10, -9) &&
            fabs(K(5,0)+(527.201615246513)) < pow(10, -9) &&

            fabs(x(0,0)-(5.74787457143478e+006)) < pow(10, -9) &&
            fabs(x(1,0)-(2.70251088466640e+006)) < pow(10, -9) &&
            fabs(x(2,0)-(3.45916012975976e+006)) < pow(10, -9) &&
            fabs(x(3,0)-(4.37955069734624e+003)) < pow(10, -9) &&
            fabs(x(4,0)+(1.94899959605007e+003)) < pow(10, -9) &&
            fabs(x(5,0)+(5.81328923582123e+003)) < pow(10, -9) &&

            fabs(P(0,0)-(95.8592545326805e+006)) < pow(10, -9) &&
            fabs(P(1,0)+(12.6624171256939e+006)) < pow(10, -9) &&
            fabs(P(2,0)-(19.4317330975035e+006)) < pow(10, -9) &&
            fabs(P(3,0)-(38.4943314161026e+003)) < pow(10, -9) &&
            fabs(P(4,0)+(670.265126224979e+000)) < pow(10, -9) &&
            fabs(P(5,0)-(10.8798071750619e+003)) < pow(10, -9) &&

            fabs(P(0,1)+(12.6624171256939e+006)) < pow(10, -9) &&
            fabs(P(1,1)-(72.2200287178262e+006)) < pow(10, -9) &&
            fabs(P(2,1)-(43.8557602801980e+006)) < pow(10, -9) &&
            fabs(P(3,1)+(624.979596635881e+000)) < pow(10, -9) &&
            fabs(P(4,1)-(25.2930978473085e+003)) < pow(10, -9) &&
            fabs(P(5,1)-(16.5443915355394e+003)) < pow(10, -9) &&

            fabs(P(0,2)-(19.4317330975035e+006)) < pow(10, -9) &&
            fabs(P(1,2)-(43.8557602801980e+006)) < pow(10, -9) &&
            fabs(P(2,2)-(35.4123169201898e+006)) < pow(10, -9) &&
            fabs(P(3,2)-(10.7791106039500e+003)) < pow(10, -9) &&
            fabs(P(4,2)-(16.4694607061714e+003)) < pow(10, -9) &&
            fabs(P(5,2)-(14.3662657116377e+003)) < pow(10, -9) &&

            fabs(P(0,3)-(38.4943314161026e+003)) < pow(10, -9) &&
            fabs(P(1,3)+(624.979596635882e+000)) < pow(10, -9) &&
            fabs(P(2,3)-(10.7791106039500e+003)) < pow(10, -9) &&
            fabs(P(3,3)-(1.00108448021159e+003)) < pow(10, -9) &&
            fabs(P(4,3)-(75.2288553456440e-003)) < pow(10, -9) &&
            fabs(P(5,3)-(3.87670231816378e+000)) < pow(10, -9) &&

            fabs(P(0,4)+(670.265126224979e+000)) < pow(10, -9) &&
            fabs(P(1,4)-(25.2930978473085e+003)) < pow(10, -9) &&
            fabs(P(2,4)-(16.4694607061714e+003)) < pow(10, -9) &&
            fabs(P(3,4)-(75.2288553456446e-003)) < pow(10, -9) &&
            fabs(P(4,4)-(996.305795884776e+000)) < pow(10, -9) &&
            fabs(P(5,4)-(5.55345184146847e+000)) < pow(10, -9) &&

            fabs(P(0,5)-(10.8798071750619e+003)) < pow(10, -9) &&
            fabs(P(1,5)-(16.5443915355394e+003)) < pow(10, -9) &&
            fabs(P(2,5)-(14.3662657116377e+003)) < pow(10, -9) &&
            fabs(P(3,5)-(3.87670231816378e+000)) < pow(10, -9) &&
            fabs(P(4,5)-(5.55345184146848e+000)) < pow(10, -9) &&
            fabs(P(5,5)-(992.788929673052e+000)) < pow(10, -9));





    return 0;
}


int Legendre_01(){

    Matrix pnm(1,1);
    Matrix dpnm(1,1);
    Legendre(3,2,1.5,pnm,dpnm);
    _assert(fabs(pnm(0,0)-1)<1e-9&&fabs(pnm(0,1)-0)<1e-9&&fabs(pnm(0,2)-0)<1e-9&&fabs(pnm(0,3)-0)<1e-9
            && fabs(pnm(1,0)-1.72771199709346)<1e-9 && fabs(pnm(1,1)-122.520427273707e-003)<1e-9 && fabs(pnm(1,2)-0)<1e-9 && fabs(pnm(1,3)-0)<1e-9
            && fabs(pnm(2,0)-2.21928488408494)<1e-9 && fabs(pnm(2,1)-273.277720516261e-003)<1e-9 && fabs(pnm(2,2)-9.68972350089721e-003)<1e-9 && fabs(pnm(2,3)-0)<1e-9
            && fabs(pnm(3,0)-2.60610986973148)<1e-9 && fabs(pnm(3,1)-455.562127741359e-003)<1e-9 && fabs(pnm(3,2)-25.5723786332915e-003)<1e-9 && fabs(pnm(3,3)-740.342454819933e-006)<1e-9
            && fabs(dpnm(0,0)-0)<1e-9 && fabs(dpnm(0,1)-0)<1e-9 && fabs(dpnm(0,2)-0)<1e-9 && fabs(dpnm(0,3)-0)<1e-9
            && fabs(dpnm(1,0)-122.520427273707e-003)<1e-9 && fabs(dpnm(1,1)+1.72771199709346)<1e-9 && fabs(dpnm(1,2)-0)<1e-9 && fabs(dpnm(1,3)-0)<1e-9
            && fabs(dpnm(2,0)-473.330896510772e-003)<1e-9 && fabs(dpnm(2,1)+3.83422445220383)<1e-9 && fabs(dpnm(2,2)+273.277720516261e-003)<1e-9 && fabs(dpnm(2,3)-0)<1e-9
            && fabs(dpnm(3,0)-1.11589475910294)<1e-9 && fabs(dpnm(3,1)+6.34320591363856)<1e-9 && fabs(dpnm(3,2)+719.400239063022e-003)<1e-9 && fabs(dpnm(3,3)+31.3196395804077e-003)<1e-9);
    return 0;
}

int Legendre_02(){

    Matrix pnm(1,1);
    Matrix dpnm(1,1);
    Legendre(2,3,6.5,pnm,dpnm);
    _assert(fabs(pnm(0,0)-1)<1e-9&&fabs(pnm(0,1)-0)<1e-9&&fabs(pnm(0,2)-0)<1e-9&&fabs(pnm(0,3)-0)<1e-9
            && fabs(pnm(1,0)-372.598749091708e-003)<1e-9 && fabs(pnm(1,1)-1.69149938580400)<1e-9 && fabs(pnm(1,2)-0)<1e-9 && fabs(pnm(1,3)-0)<1e-9
            && fabs(pnm(2,0)+962.817522589578e-003)<1e-9 && fabs(pnm(2,1)-813.649968127449e-003)<1e-9 && fabs(pnm(2,2)-1.84687740458339)<1e-9 && fabs(pnm(2,3)-0)<1e-9
            && fabs(dpnm(0,0)-0)<1e-9 && fabs(dpnm(0,1)-0)<1e-9 && fabs(dpnm(0,2)-0)<1e-9 && fabs(dpnm(0,3)-0)<1e-9
            && fabs(dpnm(1,0)-1.691499385804)<1e-9 && fabs(dpnm(1,1)+372.598749091708e-003)<1e-9 && fabs(dpnm(1,2)-0)<1e-9 && fabs(dpnm(1,3)-0)<1e-9
            && fabs(dpnm(2,0)-1.40928308437354)<1e-9 && fabs(dpnm(2,1)-3.51452627212613)<1e-9 && fabs(dpnm(2,2)+813.649968127449e-003)<1e-9 && fabs(dpnm(2,3)-0)<1e-9);
    return 0;
}

int Legendre_03(){

    Matrix pnm(1,1);
    Matrix dpnm(1,1);
    Legendre(1,2,-0.45,pnm,dpnm);
    _assert(fabs(pnm(0,0)-1)<1e-9&&fabs(pnm(0,1)-0)<1e-9&&fabs(pnm(0,2)-0)<1e-9&&fabs(pnm(0,3)-0)<1e-9
            && fabs(pnm(1,0)+753.382404621984e-003)<1e-9 && fabs(pnm(1,1)-1.55962013080301)<1e-9 && fabs(pnm(1,2)-0)<1e-9 && fabs(pnm(1,3)-0)<1e-9
            && fabs(dpnm(0,0)-0)<1e-9 && fabs(dpnm(0,1)-0)<1e-9 && fabs(dpnm(0,2)-0)<1e-9 && fabs(dpnm(0,3)-0)<1e-9
            && fabs(dpnm(1,0)-1.55962013080301)<1e-9 && fabs(dpnm(1,1)-753.382404621984e-003)<1e-9 && fabs(dpnm(1,2)-0)<1e-9 && fabs(dpnm(1,3)-0)<1e-9);
    return 0;
}

int AccelHarmonic_01(){
    double data[] = {30000.0,50000.0,90000.0};
    Matrix r = Matrix(3,1,data,3);
    Matrix E = Matrix::identity(3)*3;
    Matrix a = AccelHarmonic(r, E, 2, 2);
    _assert(fabs(a(0,0)-1.97948273625860e+003) < 1e-9
         && fabs(a(1,0)-3.27140377530750e+003) < 1e-9 &&
            fabs(a(2,0)+6.44134928737685e+003) < 1e-9);
    return 0;
}

int AccelHarmonic_02(){
    double data[] = {-36000,87000,321000};
    Matrix r = Matrix(3,1,data,3);
    Matrix E = Matrix::identity(3)*9;
    Matrix a = AccelHarmonic(r, E, 2, 2);
    _assert(fabs(a(0,0)-41.4678762157640) < 1e-9 && fabs(a(1,0)+100.227774131053) < 1e-9 &&
            fabs(a(2,0)+375.329335798179) < 1e-9);
    return 0;
}

int AccelHarmonic_03(){
    double data[] = {-389000, 107000, 901000};
    Matrix r = Matrix(3,1,data,3);
    Matrix E = Matrix::identity(3)*4;
    Matrix a = AccelHarmonic(r, E, 2, 2);
    _assert(fabs(a(0,0)-39.7499980028614) < 1e-5 && fabs(a(1,0)+10.9338223827362) < 1e-5&&
            fabs(a(2,0)+92.8621697554328) < 1e-5);
    return 0;
}

int Accel_01(){
    auxParam.Mjd_TT = 4.974611324287046e+04;
    auxParam.n = 10;
    auxParam.m = 10;
    auxParam.Mjd_UTC = 4.974611253472231e+04;
    double dataY[] = {1000000, 690000, 500000, 4, 5, 6};
    double x = 0.0;
    Matrix Y = Matrix(6,1,dataY,6);
    Matrix dY = Accel(x,Y);
    _assert(dY(0,0)==4 && dY(1,0)==5 && dY(2,0)==6 &&
            fabs(dY(3,0)+2.05530727273914e+003)<1e-7 &&
            fabs(dY(4,0)+1.96615200204762e+003)<1e-7 &&
            fabs(dY(5,0)+2.11191857745828e+003)<1e-7);
    return 0;

}

int Accel_02(){
    double dataY[] = {-2005000, 780000, 400000, -8, 9.7, 5};
    double x = 1.0;
    Matrix Y = Matrix(6,1,dataY,6);
    Matrix dY = Accel(x,Y);
    _assert(dY(0,0)==-8 && dY(1,0)==9.7 && dY(2,0)==5 &&
            fabs(dY(3,0)-77.8747455160612)<1e-7 &&
            fabs(dY(4,0)+46.5207890959991)<1e-7 &&
            fabs(dY(5,0)+12.8025487799574)<1e-7);
    return 0;

}

int Accel_03(){
    auxParam.m = 10;
    auxParam.n = 10;
    auxParam.Mjd_TT = 4.974611324287046e+04;
    double x=0.0;
    double dataY[]={6.14730428873136e+006,2.49821609757119e+006,2.87280805359544e+006,3.76462899474253e+003,-2.21784494072807e+003,-6.14147100738888e+003};
    Matrix Y(6,1,dataY,6);
    Matrix dY=Accel(x,Y);
    _assert(fabs(dY(0,0)-3.76462899474253e+003) < 1e-9 &&
            fabs(dY(1,0)+2.21784494072807e+003) < 1e-9 &&
            fabs(dY(2,0)+6.14147100738888e+003) < 1e-9 &&
            fabs(dY(3,0)+6.48314260981239) < 1e-9 &&
            fabs(dY(4,0)+2.63478439462459) < 1e-9 &&
            fabs(dY(5,0)+3.03749022533030) < 1e-9);

    return 0;
}

int G_AccelHarmonic_01(){
    double dataR[] = {5.58174858250217e+006,2.77270744874438e+006,3.67162621854941e+006};
    double dataE[] = {-976.675978197985e-003, 214.718055826991e-003 ,-436.019365249908e-006,
                      -214.718017127451e-003, -976.676074804266e-003, -134.260146426503e-006,
                      -454.677759806673e-006, -37.5074463084540e-006, 999.999895930658e-003};
    Matrix r= Matrix(3,1,dataR,3);
    Matrix E = Matrix(3,3,dataE,9);
    Matrix G = G_AccelHarmonic(r,E,10,10);
    _assert(fabs(G(0,0)-825.193994913320e-009)<1e-9 && fabs(G(0,1)-932.874952397356e-009)<1e-9 && fabs(G(0,2)-1.24046938587696e-006)<1e-9 &&
            fabs(G(1,0)-932.874954173712e-009)<1e-9 && fabs(G(1,1)+589.285963137343e-009)<1e-9 && fabs(G(1,2)-616.245046991537e-009)<1e-9 &&
            fabs(G(2,0)-1.24046938942968e-006)<1e-9 && fabs(G(2,1)-616.245047879715e-009)<1e-9 && fabs(G(2,2)+235.908035328691e-009)<1e-9);
    return 0;
}

int G_AccelHarmonic_02(){
    double dataR[] = {1.5e7, -2e7, 3e7};
    Matrix r= Matrix(3,1,dataR,3);
    Matrix E = Matrix::identity(3);
    Matrix G = G_AccelHarmonic(r,E,5,3);
    _assert(fabs(G(0,0)+3.73073678894276e-009)<1e-9 && fabs(G(0,1)+3.94918513169085e-009)<1e-9 && fabs(G(0,2)-5.92462974513541e-009)<1e-9 &&
            fabs(G(1,0)+3.94918520107979e-009)<1e-9 && fabs(G(1,1)+1.42703887595630e-009)<1e-9 && fabs(G(1,2)+7.89951515312026e-009)<1e-9 &&
            fabs(G(2,0)-5.92462984227993e-009)<1e-9 && fabs(G(2,1)+7.89951512536469e-009)<1e-9 && fabs(G(2,2)-5.15777573428799e-009)<1e-9);
    return 0;
}

int G_AccelHarmonic_03(){
    double dataR[] = {-2.5e6, -4e6, 4.3e6};
    Matrix r= Matrix(3,1,dataR,3);
    Matrix E = Matrix::identity(3)*6;
    Matrix G = G_AccelHarmonic(r,E,7,8);
    _assert(fabs(G(0,0)+137.902919639110e-009)<1e-9 && fabs(G(0,1)-188.097927455999e-009)<1e-9 && fabs(G(0,2)+202.235550084140e-009)<1e-9 &&
            fabs(G(1,0)-188.097926789865e-009)<1e-9 && fabs(G(1,1)-45.4926816217949e-009)<1e-9 && fabs(G(1,2)+323.577053595869e-009)<1e-9 &&
            fabs(G(2,0)+202.235549418006e-009)<1e-9 && fabs(G(2,1)+323.577053151780e-009)<1e-9 && fabs(G(2,2)-92.4102381283376e-009)<1e-9);
    return 0;
}

int VarEqn_01(){

    double x=0;
    auxParam.Mjd_TT=49746.11085861109313555062;
    auxParam.Mjd_UTC=49746.11015046294778585434;
    Matrix yPhi(42,1);
    yPhi(0,0)=5.58174847185176e+006;
    yPhi(1,0)=2.77270749305686e+006;
    yPhi(2,0)=3.67162635502568e+006;
    yPhi(3,0)=4.60030935493462e+003;
    yPhi(4,0)=-1.84229840480395e+003;
    yPhi(5,0)=-5.67402389148098e+003;
    yPhi(6,0)=1.00000000000000e+000;
    yPhi(7,0)=0.00000000000000e+000;
    yPhi(8,0)=0.00000000000000e+000;
    yPhi(9,0)=0.00000000000000e+000;
    yPhi(10,0)=0.00000000000000e+000;
    yPhi(11,0)=0.00000000000000e+000;
    yPhi(12,0)=0.00000000000000e+000;
    yPhi(13,0)=1.00000000000000e+000;
    yPhi(14,0)=0.00000000000000e+000;
    yPhi(15,0)=0.00000000000000e+000;
    yPhi(16,0)=0.00000000000000e+000;
    yPhi(17,0)=0.00000000000000e+000;
    yPhi(18,0)=0.00000000000000e+000;
    yPhi(19,0)=0.00000000000000e+000;
    yPhi(20,0)=1.00000000000000e+000;
    yPhi(21,0)=0.00000000000000e+000;
    yPhi(22,0)=0.00000000000000e+000;
    yPhi(23,0)=0.00000000000000e+000;
    yPhi(24,0)=0.00000000000000e+000;
    yPhi(25,0)=0.00000000000000e+000;
    yPhi(26,0)=0.00000000000000e+000;
    yPhi(27,0)=1.00000000000000e+000;
    yPhi(28,0)=0.00000000000000e+000;
    yPhi(29,0)=0.00000000000000e+000;
    yPhi(30,0)=0.00000000000000e+000;
    yPhi(31,0)=0.00000000000000e+000;
    yPhi(32,0)=0.00000000000000e+000;
    yPhi(33,0)=0.00000000000000e+000;
    yPhi(34,0)=1.00000000000000e+000;
    yPhi(35,0)=0.00000000000000e+000;
    yPhi(36,0)=0.00000000000000e+000;
    yPhi(37,0)=0.00000000000000e+000;
    yPhi(38,0)=0.00000000000000e+000;
    yPhi(39,0)=0.00000000000000e+000;
    yPhi(40,0)=0.00000000000000e+000;
    yPhi(41,0)=1.00000000000000e+000;

    Matrix yPhip=VarEqn(x,yPhi);
    //yPhip.print();
    _assert(fabs(yPhip(0,0)-(4.60030935493462e+003)) < pow(10, -9) &&
            fabs(yPhip(1,0)+(1.84229840480395e+003)) < pow(10, -9) &&
            fabs(yPhip(2,0)+(5.67402389148098e+003)) < pow(10, -9) &&
            fabs(yPhip(3,0)+(5.87599335412780e+000)) < pow(10, -9) &&
            fabs(yPhip(4,0)+(2.91894719034911e+000)) < pow(10, -9) &&
            fabs(yPhip(5,0)+(3.87497881299000e+000)) < pow(10, -9) &&
            fabs(yPhip(6,0)-(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(7,0)-(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(8,0)-(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(9,0)-(825.193922970868e-009)) < pow(10, -9) &&
            fabs(yPhip(10,0)-(932.874949732820e-009)) < pow(10, -9) &&
            fabs(yPhip(11,0)-(1.24046941163414e-006)) < pow(10, -9) &&
            fabs(yPhip(12,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(13,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(14,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(15,0)-(932.874947956464e-009)) < pow(10, -9) &&
            fabs(yPhip(16,0)+(589.285948038309e-009)) < pow(10, -9) &&
            fabs(yPhip(17,0)-(616.245080298228e-009)) < pow(10, -9) &&
            fabs(yPhip(18,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(19,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(20,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(21,0)-(1.24046940896960e-006)) < pow(10, -9) &&
            fabs(yPhip(22,0)-(616.245078521871e-009)) < pow(10, -9) &&
            fabs(yPhip(23,0)+(235.907972712113e-009)) < pow(10, -9) &&
            fabs(yPhip(24,0)-(1.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(25,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(26,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(27,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(28,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(29,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(30,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(31,0)-(1.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(32,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(33,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(34,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(35,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(36,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(37,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(38,0)-(1.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(39,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(40,0)+(0.00000000000000e+000)) < pow(10, -9) &&
            fabs(yPhip(41,0)+(0.00000000000000e+000)) < pow(10, -9));

    return 0;
}

int AzElPa_01(){
    double dataS[] = {1,2,3};
    double Az, El;
    Matrix s = Matrix(3,1,dataS,3);
    Matrix dAds = Matrix(1,1);
    Matrix dEds = Matrix(1,1);
    AzElPa(s,Az,El,dAds,dEds);
    _assert(fabs(Az-0.46365)<1e-5 && fabs(El-0.93027)< 1e-5 &&
            fabs(dAds(0,0)-400.000000000000e-003)<1e-6 && fabs(dAds(1,0)+200.000000000000e-003)<1e-6 && fabs(dAds(2,0)-0)<1e-6 &&
            fabs(dEds(0,0)+95.8314847499910e-003)<1e-6 && fabs(dEds(1,0)+191.662969499982e-003)<1e-6 && fabs(dEds(2,0)-159.719141249985e-003)<1e-6);
    return 0;
}

int AzElPa_02(){
    double dataS[] = {-6.37,10.24,34.33};
    double Az, El;
    Matrix s = Matrix(3,1,dataS,3);
    Matrix dAds = Matrix(1,1);
    Matrix dEds = Matrix(1,1);
    AzElPa(s,Az,El,dAds,dEds);
    _assert(fabs(Az-5.7267)<1e-5 && fabs(El- 1.233)< 1e-4 &&
            fabs(dAds(0,0)-70.4097033372412e-003)<1e-6 && fabs(dAds(1,0)-43.7997861580299e-003)<1e-6 && fabs(dAds(2,0)-0)<1e-6 &&
            fabs(dEds(0,0)-13.6961017789997e-003)<1e-6 && fabs(dEds(1,0)+22.0169673809980e-003)<1e-6 && fabs(dEds(2,0)-9.10859057132675e-003)<1e-6);
    return 0;
}

int AzElPa_03(){
    double dataS[] = {1.04025689965288e+006,-446.480305965849e+003,752.190746199526e+003};
    double Az, El;
    Matrix s = Matrix(3,1,dataS,3);
    Matrix dAds = Matrix(1,1);
    Matrix dEds = Matrix(1,1);
    AzElPa(s,Az,El,dAds,dEds);
    _assert(fabs(Az-1.97622068644368)<1e-5 && fabs(El-586.476946459832e-003)< 1e-5 &&
            fabs(dAds(0,0)+348.410139988782e-009)<1e-6 && fabs(dAds(1,0)+811.762685138629e-009)<1e-6 && fabs(dAds(2,0)-0)<1e-6 &&
            fabs(dEds(0,0)+374.181622245686e-009)<1e-6 && fabs(dEds(1,0)-160.599487725387e-009 )<1e-6 && fabs(dEds(2,0)-612.809350568239e-009)<1e-6);
    return 0;
}

int Sign_01(){
    _assert(sign_(-1,2) == 1);
    return 0;
}

int Sign_02(){
    _assert(sign_(-10.312,23.112) == 10.312);
    return 0;
}

int Sign_03(){
    _assert(sign_(4.3213,-4.3213) == -4.3213);
    return 0;
}

int DEInteg_01(){
    auxParam.Mjd_TT = 49746.11199287025374360382556915283203125;
    double t=0.0000000000000000000000000000000000000000;
    double tout=-134.9999919533729553222656250000000000000000;
    double relerr=0.0000000000001000000000000000030373745563;
    double abserr=0.0000009999999999999999547481118258862587;
    double n_eqn=6.000000000000000000000000000000000000000;
    Matrix y(6,1);
    y(0,0)=6147304.2887313608080148696899414062500000000000;
    y(1,0)=2498216.0975711927749216556549072265625000000000;
    y(2,0)=2872808.0535954418592154979705810546875000000000;
    y(3,0)=3764.6289947425334503350313752889633178710938;
    y(4,0)=-2217.8449407280654668284114450216293334960938;
    y(5,0)=-6141.4710073888818442355841398239135742187500;

    Matrix Y=DEInteg((*Accel),t,tout,relerr,abserr,n_eqn,y);
    Y.print();

    _assert(fabs(Y(0,0)-(5581749.4045084258541464805603027343750000000000)) < pow(10, -10) &&
            fabs(Y(1,0)-(2772707.9559171875007450580596923828125000000000)) < pow(10, -10) &&
            fabs(Y(2,0)-(3671626.9762915470637381076812744140625000000000)) < pow(10, -10) &&
            fabs(Y(3,0)-(4600.3099223832241477794013917446136474609375)) < pow(10, -10) &&
            fabs(Y(4,0)+(1842.2988221636537673475686460733413696289062)) < pow(10, -10) &&
            fabs(Y(5,0)+(5674.0250624161899395403452217578887939453125)) < pow(10, -10));

    return 0;
}

int all_tests()
{
    Cnm = new double*[361];
    for (int i=0; i<361; i++){
        Cnm[i] = new double[361];
    }
    Snm = new double*[361];
    for (int i=0; i<361; i++){
        Snm[i] = new double[361];
    }
    double temp[6];

    ifstream fid("data/egm.txt");
    if (!fid.is_open()){
        cerr << "ERROR: eegm.txt NOT FOUND" << endl;
    }

    for (int n = 0; n< 360; n++) {
        for (int j = 0; j <= n; j++) {
            for (int i = 0; i < 6; i++) {
                fid >> temp[i];
            }
            Cnm[n][j] = temp[2];
            Snm[n][j] = temp[3];
        }
    }

    fid.close();

    eopdata = new double*[19716];
    for (int i = 0; i < 19716; i++) {
        eopdata[i] = new double[13];
    }

    ifstream fid2("data/eop19620101.txt");
    if (!fid2.is_open()){
        cerr << "ERROR: eop19620101.txt NOT FOUND" << endl;
    }

    /*
     ----------------------------------------------------------------------------------------------------
     | Date MJD x y UT1-UTC LOD dPsi dEpsilon dX dY DAT
     |(0h UTC) " " s s " " " " s
     ----------------------------------------------------------------------------------------------------
    */
    for (int i = 0; i < 19716; i++) {
        for (int j = 0; j < 13; j++) {
            fid2 >> eopdata[i][j];
        }
    }

    fid2.close();

    _verify(Mjday_01);
    _verify(Mjday_02);
    _verify(Mjday_03);
    _verify(Mjday_04);
    _verify(Position_01);
    _verify(Position_02);
    _verify(Position_03);
    _verify(MOblq_01);
    _verify(MOblq_02);
    _verify(MOblq_03);
    _verify(Frac_01);
    _verify(Frac_02);
    _verify(Frac_03);
    _verify(R_x_01);
    _verify(R_x_02);
    _verify(R_x_03);
    _verify(R_x_04);
    _verify(R_y_01);
    _verify(R_y_02);
    _verify(R_y_03);
    _verify(R_y_04);
    _verify(R_z_01);
    _verify(R_z_02);
    _verify(R_z_03);
    _verify(R_z_04);
    _verify(LTC_01);
    _verify(LTC_02);
    _verify(LTC_03);
    _verify(PrecMatrix_01);
    _verify(PrecMatrix_02);
    _verify(PrecMatrix_03);
    _verify(PoleMatrix_01);
    _verify(PoleMatrix_02);
    _verify(PoleMatrix_03);
    _verify(Gmst_01);
    _verify(Gmst_02);
    _verify(Gmst_03);
    _verify(NutAngles_01);
    _verify(NutAngles_02);
    _verify(NutAngles_03);
    _verify(NutMatrix_01);
    _verify(NutMatrix_02);
    _verify(NutMatrix_03);
    _verify(EqnEquinox_01);
    _verify(EqnEquinox_02);
    _verify(EqnEquinox_03);
    _verify(Gast_01);
    _verify(Gast_02);
    _verify(Gast_03);
    _verify(Timediff_01);
    _verify(Timediff_02);
    _verify(Timediff_03);
    _verify(GHAMatrix_01);
    _verify(Cross_01);
    _verify(Cross_02);
    _verify(Cross_03);
    _verify(Norm_01);
    _verify(Norm_02);
    _verify(Norm_03);
    _verify(Dot_01);
    _verify(Dot_02);
    _verify(Dot_03);
    _verify(Transpose_01);
    _verify(Transpose_02);
    _verify(Doubler_01);
    _verify(Doubler_02);
    _verify(IERS_01);
    _verify(Geodetic_01);
    _verify(Geodetic_02);
    _verify(Geodetic_03);
    _verify(Anglesdr_01);
    _verify(TimeUpdate_01);
    _verify(TimeUpdate_02);
    _verify(TimeUpdate_03);
    _verify(MeasUpdate_01);
    _verify(MeasUpdate_02);
    _verify(MeasUpdate_03);
    _verify(Legendre_01);
    _verify(Legendre_02);
    _verify(Legendre_03);
    _verify(AccelHarmonic_01);
    _verify(AccelHarmonic_02);
    _verify(AccelHarmonic_03);
    _verify(Accel_01);
    _verify(Accel_02);
    _verify(Accel_03);
    _verify(G_AccelHarmonic_01);
    _verify(G_AccelHarmonic_02);
    _verify(G_AccelHarmonic_03);
    _verify(Sign_01);
    _verify(Sign_02);
    _verify(Sign_03);
    _verify(AzElPa_01);
    _verify(AzElPa_02);
    _verify(AzElPa_03);
    _verify(VarEqn_01);
    return 0;
}

int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    // Free memory usage
    for (int i=0; i<361; i++){
        delete[] Cnm[i];
        delete[] Snm[i];
    }
    delete[] Cnm;
    delete[] Snm;

    for (int i=0; i<19716; i++){
        delete[] eopdata[i];
    }
    delete[] eopdata;
    return result != 0;
}




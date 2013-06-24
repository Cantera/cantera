#include <fstream>
#include <cmath>
#include <sstream>

#include "cantera/base/PrintCtrl.h"
#include "base/LogPrintCtrl.h"
#include "cantera/base/global.h"

using namespace Cantera;
using namespace std;

bool dnequl(const double a, const double b, int isig)
{

    double atol = fabs(a) + fabs(b);
    double asig = pow(10.0, -isig);
    if (asig > 1.0E-6) {
        asig = 1.0E-6;
    }

    double c = fabs(a - b);
    if (c > atol * asig) {
        return true;
    }
    return false;
}

void doCheck(const double a, const double b, int isig)
{
    if (dnequl(a, b, isig)) {
        printf("error dg = %20.13e ne %g \n", a, b);
    } else {
        printf("good: dg = %20.13e eq %g \n", a, b);
    }
}


void doLogger()
{

    LogPrintCtrl ppc;
    ostream ffs(cout.rdbuf());
    double d, dg;
    int nsig, p, w, wMax, wMin;
    int Ncrop10;

    d = 0.12345;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 0.12, nsig);

    d = -1.0375E15;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.03E15, nsig);


    d = 1.0345E-270;
    nsig = 5;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0345E-270, nsig);

    d = 1.0345E-15;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.03E-15, nsig);

    d = 1.0345E-15;
    nsig = 5;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0345E-15, nsig);

    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0E-15, nsig);

    d = -1.0345E-15;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0E-15, nsig);

    d = -1.0345E-15;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.03E-15, nsig);

    d = -1.0345E-15;
    nsig = 5;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0345E-15, nsig);

    d = -1.0345E-15;
    nsig = 7;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0345E-15, nsig);


    d = 1.0345E15;
    nsig = 7;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0345E15, nsig);

    d = -1.0345E15;
    nsig = 7;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0345E15, nsig);

    d = -1.0305E15;
    nsig = 1;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0E15, nsig);


    d = -10.E-3;
    nsig = 4;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -10E-3, nsig);


    d = 1.2345;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.2, nsig);


    d = 12.345;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 12., nsig);

    d = 0.9999;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 0.99, nsig);


    d = 1.1234;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.12, nsig);

    Ncrop10 = -15;
    d = 1.1234E-15;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 1.0E-15, nsig);


    Ncrop10 = -16;
    d = 1.1234E-15;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 1.1E-15, nsig);



    Ncrop10 = -14;
    d = 1.1234E-15;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 0.0, nsig);

    Ncrop10 = 6;
    d = 1.1234E8;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 1.12E8, nsig);

    ppc.setNdec(-15);
    d = 1.1234E-14;
    p = 5;
    w = 20;
    ppc.pr_de_c10(d, p, w);
    ffs << " -- should be 1.1000E-14" << endl;


    ppc.setNdec(-14);
    d = 1.1234E-14;
    p = 5;
    w = 20;
    ppc.pr_de_c10(d, p, w);
    ffs << " -- should be 1.0000E-14" << endl;


    ppc.setNdec(-13);
    d = 1.1234E-14;
    p = 5;
    wMin = 20;
    ffs << "\"";
    ppc.pr_de_c10(d, p, wMin);
    ffs << "\" -- should be \"          0.0000E+00\" " << endl;


    ppc.setNdec(-20);
    d = 1.1234567E-14;
    p = 5;
    w = 20;
    ffs << "\"";
    ppc.pr_de_c10(d, p, w);
    ffs << "\" -- should be 1.1235E-14" << endl;


    ppc.setNdec(-20);
    d = 1.1234567E-14;
    p = 3;
    w = 17;
    ffs << "\"";
    ppc.pr_de_c10(d, p, w);
    ffs << "\" -- should be \"         1.12E-14\"" << endl;



    ppc.setNdec(-20);
    d = 1.1234567E-14;
    p = 5;
    wMin = 3;
    wMax = 8;
    ffs << "\"";
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << "\" -- should be \"1.12E-14\" " << endl;

    ppc.setNdec(-20);
    d = 1.1234567E-4;
    p = 5;
    wMin = 3;
    wMax = 8;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be 1.12E-14 " << endl;



    ppc.setNdec(-120);
    d = 1.1234567E-109;
    p = 5;
    wMin = 3;
    wMax = 8;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be 1.1E-109 " << endl;


    ppc.setNdec(-120);
    d = 1.1234567E-19;
    p = 5;
    wMin = 3;
    wMax = 9;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be 1.123E-19 " << endl;

    ppc.setNdec(-20);
    d = -1.1234567E-19;
    p = 5;
    wMin = 3;
    wMax = 10;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be -1.00E-19 " << endl;

    ppc.setNdec(-20);
    d = -1.1234567E-19;
    p = 5;
    wMin = 3;
    wMax = 11;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be -1.1000E-19 " << endl;


    ppc.setNdec(-25);
    d = -1.1274567E-19;
    p = 3;
    wMin = 3;
    wMax = 10;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be -1.13E-19 " << endl;


    ppc.setNdec(-25);
    d = -1.127401E-19;
    double dd = ppc.cropSigDigits(d, 3);
    p = 4;
    wMin = 3;
    wMax = 10;
    ppc.pr_de_c10(dd, p, wMin, wMax);
    ffs << " -- should be -1.130E-19 " << endl;

    return;
}

int main()
{
    suppress_deprecation_warnings();
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    // How to connect to a file:
    // ofstream fff("redirect.txt");
    // ostream ffs(fff.rdbuf());

    // How to connect to a string:
    // std::ostringstream os;
    // ostream ffss(os.rdbuf());

    ostream ffs(cout.rdbuf());

    PrintCtrl ppc(ffs);

    double d, dg;
    int nsig, p, w, wMax, wMin;
    int Ncrop10;

    d = 0.12345;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 0.12, nsig);

    d = -1.0375E15;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.03E15, nsig);


    d = 1.0345E-270;
    nsig = 5;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0345E-270, nsig);

    d = 1.0345E-15;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.03E-15, nsig);

    d = 1.0345E-15;
    nsig = 5;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0345E-15, nsig);

    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0E-15, nsig);

    d = -1.0345E-15;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0E-15, nsig);

    d = -1.0345E-15;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.03E-15, nsig);

    d = -1.0345E-15;
    nsig = 5;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0345E-15, nsig);

    d = -1.0345E-15;
    nsig = 7;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0345E-15, nsig);


    d = 1.0345E15;
    nsig = 7;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.0345E15, nsig);

    d = -1.0345E15;
    nsig = 7;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0345E15, nsig);


    d = -1.0305E15;
    nsig = 1;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -1.0E15, nsig);


    d = -10.E-3;
    nsig = 4;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, -10E-3, nsig);


    d = 1.2345;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.2, nsig);


    d = 12.345;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 12., nsig);

    d = 0.9999;
    nsig = 2;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 0.99, nsig);


    d = 1.1234;
    nsig = 3;
    dg = ppc.cropSigDigits(d, nsig);
    doCheck(dg, 1.12, nsig);


    Ncrop10 = -15;
    d = 1.1234E-15;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 1.0E-15, nsig);


    Ncrop10 = -16;
    d = 1.1234E-15;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 1.1E-15, nsig);



    Ncrop10 = -14;
    d = 1.1234E-15;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 0.0, nsig);

    Ncrop10 = 6;
    d = 1.1234E8;

    dg = ppc.cropAbs10(d, Ncrop10);
    doCheck(dg, 1.12E8, nsig);

    ppc.setNdec(-15);
    d = 1.1234E-14;
    p = 5;
    w = 20;
    ppc.pr_de_c10(d, p, w);
    ffs << " -- should be 1.1000E-14" << endl;


    ppc.setNdec(-14);
    d = 1.1234E-14;
    p = 5;
    w = 20;
    ppc.pr_de_c10(d, p, w);
    ffs << " -- should be 1.0000E-14" << endl;


    ppc.setNdec(-13);
    d = 1.1234E-14;
    p = 5;
    wMin = 20;
    ffs << "\"";
    ppc.pr_de_c10(d, p, wMin);
    ffs << "\" -- should be \"          0.0000E+00\" " << endl;


    ppc.setNdec(-20);
    d = 1.1234567E-14;
    p = 5;
    w = 20;
    ffs << "\"";
    ppc.pr_de_c10(d, p, w);
    ffs << "\" -- should be 1.1235E-14" << endl;


    ppc.setNdec(-20);
    d = 1.1234567E-14;
    p = 3;
    w = 17;
    ffs << "\"";
    ppc.pr_de_c10(d, p, w);
    ffs << "\" -- should be \"         1.12E-14\"" << endl;



    ppc.setNdec(-20);
    d = 1.1234567E-14;
    p = 5;
    wMin = 3;
    wMax = 8;
    ffs << "\"";
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << "\" -- should be \"1.12E-14\" " << endl;



    ppc.setNdec(-20);
    d = 1.1234567E-4;
    p = 5;
    wMin = 3;
    wMax = 8;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be 1.12E-14 " << endl;



    ppc.setNdec(-120);
    d = 1.1234567E-109;
    p = 5;
    wMin = 3;
    wMax = 8;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be 1.1E-109 " << endl;


    ppc.setNdec(-120);
    d = 1.1234567E-19;
    p = 5;
    wMin = 3;
    wMax = 9;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be 1.123E-19 " << endl;

    ppc.setNdec(-20);
    d = -1.1234567E-19;
    p = 5;
    wMin = 3;
    wMax = 10;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be -1.00E-19 " << endl;

    ppc.setNdec(-20);
    d = -1.1234567E-19;
    p = 5;
    wMin = 3;
    wMax = 11;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be -1.1000E-19 " << endl;


    ppc.setNdec(-25);
    d = -1.1274567E-19;
    p = 3;
    wMin = 3;
    wMax = 10;
    ppc.pr_de_c10(d, p, wMin, wMax);
    ffs << " -- should be -1.13E-19 " << endl;


    ppc.setNdec(-25);
    d = -1.127401E-19;
    double dd = ppc.cropSigDigits(d, 3);
    p = 4;
    wMin = 3;
    wMax = 10;
    ppc.pr_de_c10(dd, p, wMin, wMax);
    ffs << " -- should be -1.130E-19 " << endl;

    // Example of attaching to a stringstream and then outputting
    // from the string.
    std::ostringstream os;
    ostream ffss(os.rdbuf());
    PrintCtrl ppss(ffss);

    ppss.setNdec(-25);
    d = -1.127401E-19;
    dd = ppss.cropSigDigits(d, 3);
    p = 4;
    wMin = 3;
    wMax = 10;
    ppss.pr_de_c10(dd, p, wMin, wMax);
    ffs << "\"";
    ffs << os.str();
    ffs << "\" -- should be \"-1.130E-19\" " << endl;

    doLogger();

}

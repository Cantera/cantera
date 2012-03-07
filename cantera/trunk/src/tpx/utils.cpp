#include "subs.h"
#include "cantera/tpx/utils.h"

using std::string;

namespace tpx
{

static string lowercase(string s)
{
    size_t i, n = s.size();
    string lc(s);
    for (i = 0; i < n; i++) {
        lc[i] = (char) tolower(s[i]);
    }
    return lc;
}

Substance* GetSubstanceByName(string name)
{
    string lcname = lowercase(name);
    if (lcname == "water") {
        return new water;
    } else if (lcname == "nitrogen") {
        return new nitrogen;
    } else if (lcname == "methane") {
        return new methane;
    } else if (lcname == "hydrogen") {
        return new hydrogen;
    } else if (lcname == "oxygen") {
        return new oxygen;
    } else if (lcname == "hfc134a") {
        return new HFC134a;
    } else if (lcname == "rk") {
        return new RedlichKwong;
    } else if (lcname == "carbondioxide") {
        return new CarbonDioxide;
    } else if (lcname == "heptane") {
        return new Heptane;
    } else {
        return 0;
    }
}

Substance* GetSub(int isub)
{
    if (isub == 0) {
        return new water;
    } else if (isub == 1) {
        return new nitrogen;
    } else if (isub == 2) {
        return new methane;
    } else if (isub == 3) {
        return new hydrogen;
    } else if (isub == 4) {
        return new oxygen;
    } else if (isub == 5) {
        return new HFC134a;
    } else if (isub == 6) {
        return new RedlichKwong;
    } else if (isub == 7) {
        return new CarbonDioxide;
    } else if (isub == 8) {
        return new Heptane;
    } else {
        return 0;
    }
}

}


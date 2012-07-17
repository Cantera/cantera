#ifndef CT_CK2CT_H
#define CT_CK2CT_H

#include <iostream>
#include <string>
#include <cstdlib>

//#include "cantera/base/ctml.h"

namespace ckr
{
class CKReader;
}

namespace pip
{

void ck2ct(std::string idtag, ckr::CKReader& r);

int convert_ck(const char* in_file, const char* db_file,
               const char* tr_file, const char* id_tag, bool debug, bool validate);

}

#endif


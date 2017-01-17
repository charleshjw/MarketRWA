// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

//#define BOOST_DATE_TIME_NO_LIB // because of the Linker problem "cannot open file libboost_date_time-vc ... Before include headers
#include <stdio.h>
#include <sstream>
#include <tchar.h>
//#include <algorithm>

#include <map>
#include <boost/assign.hpp>
#include <boost/algorithm/string/trim.hpp>

using namespace boost::algorithm;
using namespace boost::assign;

#include "Currency.h"
#include "portfolio.h"
#include "Factor.h"
#include "CoreFactor.h"
#include "ParamControl.h"
#include "ErrorLog.h"
#include "MathRoutines.h"
#include "utilities.h"
#include "CCAR_info.h"
#include "CCAR_shift.h"
#include "ShiftScenInfo.h"
std::vector<boost::shared_ptr<Factor> >::iterator findFactor(FactType iType, std::string subType, std::string cur, std::vector<boost::shared_ptr<Factor> >& facts);

// TODO: reference additional headers your program requires here

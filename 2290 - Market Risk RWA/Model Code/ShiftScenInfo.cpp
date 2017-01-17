#include "stdafx.h"

ShiftScenInfo::ShiftScenInfo(std::string nm1, std::string nm2, unsigned int n1, unsigned int n2, ShiftType sh, double v) : _s_name1(nm1), _s_name2(nm2), _tenor1(n1), _tenor2(n2), _shiftType(sh), _shift(v) {
	std::vector<std::string> values;

	boost::algorithm::split(values, nm1, is_any_of(","));
	if(values[0] == "factorGB") {
		_s_name1 = values[0].c_str();
		_s_name2 = _s_name1;
		_tenor1 = atoi(values[1].c_str());
		boost::algorithm::split(values, nm2, is_any_of(","));
		_tenor2 = atoi(values[1].c_str());
	}
}

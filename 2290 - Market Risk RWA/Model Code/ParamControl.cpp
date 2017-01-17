#include "stdafx.h"

using namespace std;

void ParamControl::readFromFile(std::string fnm) {
	ifstream inpf;
	inpf.open(fnm.c_str());
	string buf("");
	std::vector<string> values;

	while(getline(inpf, buf)) {
		boost::algorithm::split(values, buf, boost::is_any_of(","));
		if(buf.find("minDatesCount")!=std::string::npos) 
			_minDatesCount = atoi(values[1].c_str());

		if(buf.find("num_pca_ir_curve")!=std::string::npos) 
			_num_pca_ir_curve = atoi(values[1].c_str());
		if(buf.find("num_pca_fx_forward")!=std::string::npos) 
			_num_pca_fx_forward = atoi(values[1].c_str());
		if(buf.find("num_pca_eqt_spot")!=std::string::npos) 
			_num_pca_eqt_spot = atoi(values[1].c_str());
		if(buf.find("num_pca_eqt_vol")!=std::string::npos) 
			_num_pca_eqt_vol = atoi(values[1].c_str());
		if(buf.find("num_pca_fx_spot")!=std::string::npos) 
			_num_pca_fx_spot = atoi(values[1].c_str());
		if(buf.find("num_pca_fx_vol")!=std::string::npos) 
			_num_pca_fx_vol = atoi(values[1].c_str());

		if(buf.find("minNumRetsForRegr")!=std::string::npos) 
			_minNumRetsForRegr = atoi(values[1].c_str());
		if(buf.find("maxAllowedGap")!=std::string::npos) 
			_maxAllowedGap = atoi(values[1].c_str());
		if(buf.find("maxDatesUnch")!=std::string::npos) 
			_maxDatesUnch = atoi(values[1].c_str());
		if(buf.find("defaultIdioVol")!=std::string::npos) 
			_defaultIdioVol = atof(values[1].c_str());
		if(buf.find("defaultIdioEquityVol")!=std::string::npos) 
			_defaultIdioEquityVol = atof(values[1].c_str());
		if(buf.find("staleDataDailyVolMin")!=std::string::npos) 
			_staleDataDailyVolMin = atof(values[1].c_str());
		if(buf.find("MC_sumulations_num")!=std::string::npos) 
			_MC_num = atoi(values[1].c_str());
		if(buf.find("decay_factor")!=std::string::npos) 
			_decay_factor = atof(values[1].c_str());
		if(buf.find("percentile")!=std::string::npos) 
			_percentile = atof(values[1].c_str());
		if(buf.find("randomness_control")!=std::string::npos) 
			_randomness_control = atoi(values[1].c_str());
		if(buf.find("largest_return")!=std::string::npos) 
			_largestAllowedReturn = atof(values[1].c_str());
		if(buf.find("shift")!=std::string::npos) 
			_shift = atof(values[1].c_str());
		if(buf.find("sim_horizon")!=std::string::npos) 
			_sim_horizon = atoi(values[1].c_str());
		if(buf.find("MC_num")!=std::string::npos) 
			_MC_num = atoi(values[1].c_str());
	}                                                                                                                     
		
	inpf.close();
}
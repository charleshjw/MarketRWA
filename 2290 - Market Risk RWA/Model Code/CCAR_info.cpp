#include "stdafx.h"

std::map<FactType, std::string> default_eod_map = map_list_of (FX_Forward, "USDInterbank") (IR_Curve, "USDInterbank") (Eqt_Spot,  "C00000117_Index_SPALNSpot") (EQ_Vol, "SP500_C00000117_EQTV")
	(IRCap_Vol, "IR02Vol") (FX_Pair_Vol, "V44Vol");
std::map<FactType, std::string>  default_vols_map = map_list_of (FX_Forward, "USDInterbank") (IR_Curve, "USDInterbank") (Eqt_Spot,  "C00000117_Index_SPALNSpot") (EQ_Vol, "SP500_C00000117_EQTV")
	(IRCap_Vol, "IR02Vol") (FX_Pair_Vol, "V44Vol");


bool CCAR_info::read_info(std::string dir) {
	bool core_rf_info_b(false), reg_2_core_cur_b(false), rfTypeCur_coreRf_b(false), core_rf_projections_b(false), readOk(false);
	std::string core_rf_info_fn(dir + "ccar_core_rf_info.csv"), reg_2_core_cur_fn(dir + "ccar_reg_2_coreCur.csv"), rfType_2_coreRf_fn(dir + "ccar_rfType_2_coreRf.csv"), core_rf_proj_fn(dir + "ccar_core_rf_projections.csv");
	//std::string core_fact_proj_fn("core_fact_projections"), core_rf_info_fn("core_ccar_rf_info.csv"), reg_2_core_cur_fn("reg_2_core_cur.csv"), factType_2_coreFact_fn("factType_2_coreFact.csv");

	_num_projections = 0;

	//core_rf_info_b =		read_core_rf_info(core_rf_info_fn);
	reg_2_core_cur_b =		read_reg_2_core_cur_map(reg_2_core_cur_fn);
	rfTypeCur_coreRf_b =	read_rfTypeCur_coreRf_map(rfType_2_coreRf_fn);
	core_rf_projections_b = read_core_rf_projections(core_rf_proj_fn);
	
	readOk = reg_2_core_cur_b && rfTypeCur_coreRf_b && core_rf_projections_b;
	return readOk;
}
bool CCAR_info::read_reg_2_core_cur_map(std::string fn) {
	std::ifstream ff;
	ff.open(fn);
	std::string buf("");
	std::vector<std::string> values;
	
	if(!ff.is_open()) 
		return false;

	while(getline(ff, buf)) {
		if(buf.find("Core Currency for CCAR") != std::string::npos) // header line
			continue;

		boost::algorithm::split(values, buf, is_any_of(",")); 
		if(values.size() >= 2)
			_reg_2_cur_mp.insert(std::make_pair<std::string, std::string>(values[0].c_str(), values[1].c_str()));
	}
	ff.close();
	return true;
}
bool CCAR_info::read_rfTypeCur_coreRf_map(std::string fn) {
	std::ifstream ff;
	ff.open(fn);
	std::string buf("");
	std::vector<std::string> values;
	std::pair<std::string, std::string> rfType_Cur_pair;
	// std::map<std::pair<std::string, std::string>, std::string> _rfTypeCur_coreRf_mp
	if(!ff.is_open()) 
		return false;

	while(getline(ff, buf)) {
		if(buf.find("Mapped to (for CCAR)") != std::string::npos) // header line
			continue;

		boost::algorithm::split(values, buf, is_any_of(",")); 
		if(values.size() >= 3) {
			rfType_Cur_pair = std::make_pair<std::string, std::string>(values[0].c_str(), values[1].c_str());
			_rfTypeCur_coreRf_mp.insert(std::make_pair<std::pair<std::string, std::string>, std::string>(rfType_Cur_pair, values[2].c_str()));
		}
	}
	ff.close();
	return true;
}

bool CCAR_info::read_core_rf_info(std::string fn) {
	
	std::ifstream ff;
	ff.open(fn);
	std::string buf("");
	std::vector<std::string> values;
	bool flag;

	if(!ff.is_open()) 
		return false;

	while(getline(ff, buf)) {
		if(buf.find("CCAR core risk factor") != std::string::npos)
			continue;

		boost::algorithm::split(values, buf, is_any_of(","));
		if(values.size()  >= 2){
			
			if(values[1].c_str() == std::string("TRUE"))
				flag = true;
			else
				flag = false;
			
			_core_fact_info.insert(std::make_pair<std::string, bool>(values[0].c_str(), flag));
		}
	}
	ff.close();
	return true;
}

bool CCAR_info::read_core_rf_projections(std::string fn) {
	std::string buf(""), prev_name(""), name;
	std::vector<std::string> values;
	int tenor;
	unsigned int i;
	std::map<unsigned int, std::vector<double> > levels;
	std::vector<double> proj;
	bool flag;

	std::ifstream ff;
	ff.open(fn);
	if(!ff.is_open()) 
		return false;

	while(getline(ff, buf)) {

		if(buf.find("Core RF") != std::string::npos) // header line
			continue;

		boost::algorithm::split(values, buf, is_any_of(","));
		if(values.size() < 4) 
			continue;
		
		name = values[0].c_str();
		if(name != prev_name && prev_name != "") {
			_eod_proj.insert(std::make_pair<std::string, std::map<unsigned int, std::vector<double> > >(prev_name, levels));
			_core_fact_info.insert(std::make_pair<std::string, bool>(prev_name, flag));
			levels.clear();
		}
		prev_name = name;

		if(_num_projections == 0)
			_num_projections = values.size() - 3;

		if(values[1].c_str() == "TRUE")
			flag = true;

		tenor = atoi(values[2].c_str());
		if(values.size() != _num_projections + 3 || tenor < 0)
			return false;

		proj.clear();
		for(i = 0; i < _num_projections; ++i)
			proj.push_back(atof(values[i + 3].c_str()));

		levels.insert(std::make_pair<unsigned int, std::vector<double> >(tenor, proj));
	}
	_eod_proj.insert(std::make_pair<std::string, std::map<unsigned int, std::vector<double> > >(prev_name, levels));
	_core_fact_info.insert(std::make_pair<std::string, bool>(prev_name, flag));

	ff.close();
	return true;
}
bool CCAR_info::find_projection(std::string name,  bool isEOD, unsigned int i_proj, std::map<double, double>& shifts) const {
	// map: (tenor, scaler)
	bool found(false);
	std::map<unsigned int, std::vector<double> > mp;
	std::map<unsigned int, std::vector<double> >::iterator itm;
	//std::map<std::string, std::pair<double, std::vector<double> > >::iterator itmap;
	// std::map<std::string, std::map<unsigned int, std::vector<double> > > _eod_map;
	std::map<std::string, std::map<unsigned int, std::vector<double> > >::const_iterator itmap;
	std::vector<double> vec;
	double tenor;
	matrix<double> mat;
	
	shifts.clear();
	if(isEOD) {
		itmap = _eod_proj.find(name);
		if(itmap != _eod_proj.end()) 
			found = true;
	}
	else {
		itmap = _vol_proj.find(name);
		if(itmap != _vol_proj.end()) 
			found = true;
	}
	if(found) {
		mp = itmap->second;
		
		for(itm = mp.begin(); itm != mp.end(); ++itm) {

			tenor = (*itm).first;
			vec = (*itm).second;
			shifts.insert(std::make_pair<double, double>(tenor, vec[i_proj])); 
		}
	}
	
	return found;
}
void CCAR_info::calc_crf_shifts(const FactSet& facs, std::vector<boost::shared_ptr<CCAR_shift> >& ccar_shifts) {
	std::string s_name, core_vol_rf(":SP500_C00000117_EQTV");

	std::map<std::string, std::map<unsigned int, std::vector<double> > >::iterator itmap;
	std::map<unsigned int, std::vector<double> > shmp;
	std::map<double, double> mp;
	std::map<unsigned int, std::vector<double> >::iterator itm;
	std::map<std::string, bool>::iterator itmap_ismult;
	std::vector<double> tenors, vals, eod;
	std::vector<unsigned int> terms;
	bool isMult(true), isLevel(true);
	matrix<double> proj, shifts;
	unsigned int j, ii, nterms;
	
	std::vector<boost::shared_ptr<Factor> >::const_iterator jtfac;
	std::vector<boost::shared_ptr<CCAR_shift> >::iterator it_shifts;

	for(itmap = _eod_proj.begin(); itmap != _eod_proj.end(); ++itmap) {
		s_name = itmap->first;
		jtfac = facs.findFactor(s_name, false);
		if(jtfac == facs.getFactors().end())
			continue;

		boost::shared_ptr<CCAR_shift> ptr(new CCAR_shift);
		ptr->setIsLevel(isLevel);
		ptr->setName(s_name);

		itmap_ismult = _core_fact_info.find(s_name);
		if(itmap_ismult != _core_fact_info.end()) {
			isMult = itmap_ismult->second;
			ptr->setIsMult(isMult);
		}
		shmp = itmap->second;
		tenors.resize(shmp.size());
		proj.resize(shmp.size(), _num_projections);
		j = 0;
		for(itm = shmp.begin(); itm != shmp.end(); ++itm) {
			tenors[j] = itm->first;
			vals = itm->second;
			for(ii = 0; ii < _num_projections; ++ii)
				proj(j, ii) = vals[ii];
			++j;
		}
		nterms = (*jtfac)->getNterms();
		terms = (*jtfac)->getTerms();
		eod = (*jtfac)->getEODVal();
		mp.clear();
		for(j = 0; j < terms.size(); ++j) 
			mp.insert(std::make_pair<double, double> (terms[j], eod[j]));
		
		vals.resize(tenors.size());
		shifts.resize(tenors.size(), _num_projections);
		for(ii = 0; ii < _num_projections; ++ii) {
			
			interpolate_linear(mp, tenors, vals);
			for(j = 0; j < tenors.size(); ++j)
				if(isMult) 
					shifts(j, ii) = proj(j, ii) / vals[j];
				else
					shifts(j, ii) = proj(j, ii) - vals[j];
		}
		ptr->setTerms(tenors);
		ptr->setShifts(shifts);
		ptr->setIsLevel(true);
		ccar_shifts.push_back(std::move(ptr));
		if(s_name == core_vol_rf) { // Vol factor
			boost::shared_ptr<CCAR_shift> ptr(new CCAR_shift);
			ptr->setName(s_name);
			ptr->setTerms(tenors);
			ptr->setShifts(shifts);
			ptr->setIsLevel(false);
			ptr->setIsMult(true);
			ccar_shifts.push_back(std::move(ptr));
		}
	}
}
/*
for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {

		itmap = core_rf_info_mp.find((*itfac)->getScenName());
		itmap_proj = eod_proj_mp.find((*itfac)->getScenName());
		if(itmap == core_rf_info_mp.end() || itmap_proj == eod_proj_mp.end()) 
			continue;
		
		(*itfac)->setCoreRfFlag(true);
		(*itfac)->setCoreRfType(itmap->second);
		eod_proj = itmap_proj->second;
		proj_ori.resize(eod_proj.size(), n_proj);
		proj_tenors.resize(eod_proj.size());
		ii = 0;
		for(itm = eod_proj.begin(); itm != eod_proj.end(); ++itm) {
			proj_tenors[ii] = itm->first;
			
			for(j = 0; j < n_proj; ++j)
				proj_ori(ii, j) = itm->second[j];
			++ii;
		}
		nterms = (*itfac)->getNterms();
		tenors.clear();
		val.resize(nterms);
		ccar_shifts.resize(nterms, n_proj);
		for(j = 0; j < nterms; ++j)
			tenors.push_back((*itfac)->getTerms()[j]);

		for(ii = 0; ii < n_proj; ++ii) {
			term_val_mp.clear();
			for(j = 0; j < eod_proj.size(); ++j) 
				term_val_mp.insert(std::make_pair<double, double> (proj_tenors[j], proj_ori(j, ii)));
			
			interpolate_linear(term_val_mp, tenors, val);
			for(j = 0; j < nterms; ++j)
				ccar_shifts(j, ii) = val[j] / (*itfac)->getEODVal()[j];
		}
		(*itfac)->setCCARshifts(ccar_shifts);
	}
	*/
/*
bool CCAR_info::find_mapping(std::string name, bool isEOD, std::string& mapped_name, std::map<unsigned int, std::vector<double> >& shift, bool& mult) {
	bool found(false);
	std::map<std::string, std::map<unsigned int, std::vector<double> > >::iterator itmap;
	std::map<std::string, std::map<unsigned int, std::vector<double> > > mp;

	if(isEOD)
		mp = _eod_map;
	else
		mp = _vol_map;

	itmap = mp.find(name);
	if(itmap != mp.end()) {
		mapped_name = itmap->first;
		shift = itmap->second;
		found = true;
	}
	return found;
}
*/
/*
bool CCAR_info::read_core_rf_projections(std::string fn) {
	std::ifstream ff;
	ff.open(fn);
	std::string buf("");
	std::vector<std::string> values;

	if(!ff.is_open()) 
		return false;

	while(getline(ff, buf)) {
		if(buf.find("Mapped to (for CCAR)")) // header line
			continue;

		boost::algorithm::split(values, buf, is_any_of(",")); 
		if(values.size() >= 3) {

		}
	}
	return true;
}
*/
/*
bool CCAR_info::read_map(std::string fn, std::map<std::string, std::map<unsigned int, std::vector<double> > >& map) {
	std::ifstream map_f;
	std::string buf("");
	std::vector<std::string> values;
	map_f.open(fn);
	if(!map_f.is_open()) 
		return false;

	while(getline(map_f, buf)) {
		boost::algorithm::split(values, buf, is_any_of(","));
		if(values.size()  >= 3){
			std::pair<std::string, double> p(values[1].c_str(), atof(values[2].c_str()));
			map.insert(std::make_pair<std::string, std::pair<std::string, double> >(values[0].c_str(), p));
		}
	}
	
	map_f.close();
	return true;
}
*/
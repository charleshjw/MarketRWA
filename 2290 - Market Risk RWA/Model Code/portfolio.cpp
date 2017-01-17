#include "stdafx.h"

extern std::ofstream eig_info, log_ff; // temp Tanya
extern std::vector<std::string> factTypeNames;

FactSet::FactSet(Factor& ff) : _systematicVar(0.0), _idioVar(0.0),_totalVar(0.0) {
	boost::shared_ptr<Factor> ptr(&ff);
	_pos.push_back(ptr);
}
FactSet::FactSet(const std::vector<boost::shared_ptr<Factor> >& facs) {
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac;
	for(itfac = facs.begin(); itfac != facs.end(); ++itfac) 
		_pos.push_back(std::move(*itfac));
}
void FactSet::readFromScenValFile(std::string fn, std::string info_fn, ValFileType type, std::vector<boost::shared_ptr<Factor> >& facs) {
	// Read valuation results from the file
	std::ifstream val_file, info_file;
	val_file.open(fn);
	info_file.open(info_fn.c_str());

	if(!val_file.is_open() || !info_file.is_open())
		return;

	std::string buf(""), name;
	unsigned int nscen, nport, i, iscen;
	double val;
	std::vector<std::string> values;
	std::vector<std::string> pfl_names;
	matrix<double> prc;

	// Read from valuation file
	// First line: Number of Scenarios, Number of FactSets
	// Port1 0 value0
	// Port1 1 value1
	// ...
	// Port2 0 value0
	// Port2 1 value1
	getline(val_file, buf);
	boost::algorithm::split(values, buf, boost::is_any_of(" ")); 
	if(values.size() != 3) 
		return;

	nscen = atoi(values[1].c_str());
	nport = atoi(values[2].c_str());
	prc.resize(nscen, nport);

	for(i = 0; i < nport; ++i) {
		for(iscen = 0; iscen < nscen; ++iscen) {
			if(!getline(val_file, buf))
				return;
			boost::algorithm::split(values, buf, boost::is_any_of(" ")); 
			if(values.size() != 3) 
				return;
			name = values[0].c_str();
			val = atof(values[2].c_str());
			if(iscen == 0)
				pfl_names.push_back(name);
			prc(iscen, i) = val;
		}
	}

	// Read from factor info file
	getline(info_file, buf); // comment, shift, shift size
	getline(info_file, buf); // header ( "scenario,factor name, tenor")
	while(getline(info_file, buf)) {
		boost::algorithm::split(values, buf, boost::is_any_of(",")); 
		if(values.size() == 4) { // Delta
			// scenname, fact hist name, tenor, val


		}
		if(values.size() == 6) { // Cross-Gammas


		}

	}
}
void FactSet::setSensitivities(std::vector<boost::shared_ptr<Factor> >& facs) {
	// Set delta to 1, include all factors to the portfolio
	unsigned int nterms;
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::vector<double> delta;

	_pos.clear();
	for(itfac= facs.begin(); itfac != facs.end(); ++itfac) {
		nterms = (*itfac)->getNterms();
		if(nterms > 0) {
			delta.resize(nterms, 1.0);
			(*itfac)->setDelta(delta);
		}
		_pos.push_back(*itfac);
	}
}
void FactSet::setSensitivities(std::string sens_fn, std::vector<boost::shared_ptr<Factor> >& facs) {

	// Read sensitivities from file, only delta so far
	std::string buf(""), name;
	unsigned int iterm;
	std::ifstream sens_file;
	std::vector<double> sens;
	std::vector<std::string> values;

	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac;
	double val(0.0);

	sens_file.open(sens_fn.c_str());
	getline(sens_file, buf);
	_pos.clear();
	if(sens_file.is_open()) {
		while(getline(sens_file, buf)) {
			boost::algorithm::split(values, buf, boost::is_any_of(",")); 
			if(values.size() == 4) {
				name = values[0];
				itfac = matchNameWithFactor(name, facs, true);
				if(itfac != facs.end()) {
					val = atof(values[2].c_str());
					iterm = atoi(values[1].c_str());
					
					(*itfac)->setDelta(iterm, val);
					//(*itfac)->setGamma(sens);
					_pos.push_back(*itfac);
				}
			}
		}
	} 
	sens_file.close();

}
double FactSet::calc_MC_Var(double percentile) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	double this_factor_var(0);;
	unsigned int nsim(_pos[0]->_sim_rets.size1()), i;
	int m;
	std::vector<double> pfl_pnl(nsim, 0), fact_pnl(nsim, 0);
	bool isPnLCalced(false);
	m = (nsim * percentile) - 1;
	if(m < 0)
		return 0.0;

	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		isPnLCalced = (*itfac)->calcPnL_MC(fact_pnl, m);
		(*itfac)->setVar(-12345.0);
		if(isPnLCalced) {
			for(i = 0; i < nsim; ++i) {
				pfl_pnl[i] += fact_pnl[i];
			}
			sort(fact_pnl.begin(), fact_pnl.end());
			this_factor_var = fact_pnl[m];
			(*itfac)->setVar(this_factor_var);
		}
	}
	
	sort(pfl_pnl.begin(), pfl_pnl.end());
	if(m < pfl_pnl.size())
		_totalVar = pfl_pnl[m];
	return _totalVar;
}
void FactSet::print_vars(std::string fn) {
	unsigned int i(0);
	std::ofstream var_out;
	std::stringstream out_file;

	var_out.open(fn.c_str());
	var_out << "All factors VAR:," << _totalVar << std::endl;
	var_out << "Term, Name, Mapping Type, Asset Class, Currency, region, EOD value, Hist Vol, Ratio of Idiosyncr and Sstematic variances, Delta, Gamma,VaR by term, VaR" << std::endl;
	for(i = 0; i != _pos.size(); ++i)
			_pos[i]->print_vars(out_file);

	var_out << out_file.rdbuf();
	var_out.close();
}

std::vector<boost::shared_ptr<Factor> >::const_iterator FactSet::findFactor(std::string name, bool isHistName) const {
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac(_pos.begin());

	if(isHistName)
		while(itfac != _pos.end() && name != (*itfac)->getHistName())
			++itfac;
	else
		while(itfac != _pos.end() && name != (*itfac)->getScenName())
			++itfac;
	return itfac;
}
void FactSet::print_sensitivities(std::string fn) {
	std::ofstream sens_file;
	std::stringstream out_file;
	sens_file.open(fn.c_str());

	if(!sens_file.is_open())
		return;
	
	unsigned int i;
	for(i = 0; i < _pos.size(); ++i)
		if(!_pos[i]->getExcludeFromScenarios())
			_pos[i]->print_sensitivities(out_file);

	sens_file << out_file.rdbuf();
	sens_file.close();
}
void FactSet::clear_CCAR_info(const CCAR_info& info, unsigned int ii) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		(*itfac)->reset_eod();
		(*itfac)->reset_rets();
	}
}
void FactSet::save_rets_eod() {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		(*itfac)->save_rets();
		(*itfac)->save_eod();
	}
}
void FactSet::reset_rets_eod() {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		(*itfac)->reset_rets();
		(*itfac)->reset_eod();
	}
}
void FactSet::output_sim_diff(std::string fn) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::ofstream ff;
	ff.open(fn);
	if(ff.is_open()) {
		ff << "h_name,s_name, nTerms,Currency,Region, Mapping Type,Factor Type,Is Lognormal,Tenor #,Tenor, Historical returns, Dist Type, Hist Vol, EOD,Stdev of Diff, sim," << std::endl;
		for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
			(*itfac)->print_sim_diff(ff);
		}
		ff.close();
	}
}
void FactSet::output_hist_rets(std::string fn) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::ofstream ff;
	ff.open(fn);
	if(ff.is_open()) {
		for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac)
			(*itfac)->print_rets(ff);
		ff.close();
	}
}
void FactSet::readEOD(std::string dir) {
	std::ifstream ff;
	std::string fn[]={"ir_zerocurves.csv", "fx_zerocurves.csv", "other_zerocurves.csv", "creditbond_zerocurves.csv", "IndexCurves.csv", "MBS_EOD_Values.csv"};  //Temp Tanya
	std::string str;
	static const std::vector<std::string> file_names(fn, fn + sizeof(fn)/sizeof(fn[0]));
	std::vector<unsigned int> ncolumns, name_column, terms;
	unsigned int i, numFiles(file_names.size());
	std::string buf(""), s_name("");
	std::vector<std::string> values;
	bool newFactorFound(false), factorFound(false);
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac(_pos.end()), jtfac;
	std::vector<double> arg, vals;
	std::map<double, double> mp;
	std::map<double, double>::iterator it, jt;
	double vv(0.0);

	ncolumns.resize(numFiles);
	name_column.resize(numFiles);

	// Temp Tanya this should be changed
	for(i = 0; i < numFiles; ++i) {
		if(file_names[i] == "ir_zerocurves.csv") {
			ncolumns[i] = 13; // is_zerocurve.csv
			name_column[i] = 2;
		}
		else {
			ncolumns[i] = 12; // fx_zerocurves.csv, other_zerocurves.csv, creditbond_zerocurves.csv, IndexCurves.csv, MBS_EOD_Values.csv
			name_column[i] = 3;
		}
	}
	for(i = 0; i < numFiles; ++i) {
		ff.open(dir + file_names[i]);
		
		while(getline(ff, buf)) {
			if(buf.find("RiskMetrics Link") != std::string::npos) {
			//if(buf.find(",,,,,,,") == std::string::npos) {
				//if((*itfac)->getHistName() == "USDMUCCCSP" || (*itfac)->getHistName() == "EURMM")
				//	int tt = 0;
				if(itfac != _pos.end())  // Process the factor
					(*itfac)->assignEODVal(mp, *this);
				mp.clear();
				newFactorFound = true;
				itfac = _pos.end();
				continue;
			}
			boost::algorithm::split(values, buf, is_any_of(",")); 
			if(newFactorFound) {
				if(values.size() > 4) {
					s_name = ":" + values[name_column[i]];
					itfac = matchNameWithFactor(s_name, _pos, false);
					newFactorFound = false;
				}
				continue;
			}
			if(itfac != _pos.end()) {
				if(values.size() == ncolumns[i]) 
					mp.insert(std::make_pair(atoi(values[ncolumns[i] - 2].c_str()), atof(values[ncolumns[i] - 1].c_str())));
			}
		}
		ff.close();
	}
	if(itfac != _pos.end())  // Process the last factor
		(*itfac)->assignEODVal(mp, *this);
}
void FactSet::readEOD_equity_fx(std::string dir) {
	std::ifstream ff;
	std::string fn[]={"murexEquityRate.csv", "cmiEquityRate.csv", "fx_spot.csv"};  //Temp Tanya
	std::string str, keyWord("");
	static const std::vector<std::string> file_names(fn, fn + sizeof(fn)/sizeof(fn[0]));
	std::vector<unsigned int> ncolumns, name_column, terms;
	unsigned int i, numFiles(file_names.size());
	std::string buf(""), s_name("");
	std::vector<std::string> values;
	bool newFactorFound(false), factorFound(false);
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac(_pos.end()), jtfac;
	std::vector<double> arg, vals;
	std::map<double, double> mp;
	std::map<double, double>::iterator it, jt;
	double vv(0.0);

	ncolumns.resize(numFiles);
	name_column.resize(numFiles);

	for(i = 0; i < numFiles; ++i) {
		ncolumns[i] = 13; 
		name_column[i] = 3;
		if(file_names[i] == "fx_spot.csv") {
			ncolumns[i] = 14; 
			name_column[i] = 3;
		}
	}
	for(i = 0; i < numFiles; ++i) {
		ff.open(dir + file_names[i]);

		keyWord = "Index";
		if(file_names[i] == "fx_spot.csv")
			keyWord = "Exchange Rate";

		while(getline(ff, buf)) {
			boost::algorithm::split(values, buf, is_any_of(","));
			if(buf.find(keyWord) != std::string::npos) {
	
				if(itfac != _pos.end())  // Process the factor
					(*itfac)->assignEODVal(mp, *this);
				mp.clear();
				//newFactorFound = true;
				str = ":" + values[name_column[i]];
				if(file_names[i] != "fx_spot.csv") {
					str = str.substr(0, str.length() - 1);
					str = str + "Return";
				}
				s_name = str;
				itfac = matchNameWithFactor(s_name, _pos, false);
				continue;
			}
			if(itfac != _pos.end() && values.size() == ncolumns[i]) {
				if(file_names[i] != "fx_spot.csv")
					mp.insert(std::make_pair(atoi(values[ncolumns[i] - 3].c_str()), atof(values[ncolumns[i] - 2].c_str())));
				else
					mp.insert(std::make_pair(atoi(values[ncolumns[i] - 2].c_str()), atof(values[ncolumns[i] - 1].c_str())));
			}
		}
		ff.close();
	}
	if(itfac != _pos.end())  // Process the last factor
		(*itfac)->assignEODVal(mp, *this);
}
void FactSet::getFactorsByTypeAndDataAval(HistDataAvailabilityType dataAvalType, FactType factType, std::vector<boost::shared_ptr<Factor> >& f_vec) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		if((*itfac)->getType() == factType && (*itfac)->getDataAvalType() == dataAvalType && (*itfac)->getNterms() > 0)
			f_vec.push_back(*itfac);
	}
			
}

void FactSet::map_2_regress(ParamControl& params, std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Spot, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Spot, 
	std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_IR_Curve,
	std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Forward, 	std::vector<date>& dts, std::map<std::string, std::vector<std::string> >& reg_2_cur_map, std::vector<date>& dts_common) {
	bool mapped;
	
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac=_pos.begin(); itfac != _pos.end(); ++itfac) {
		if(*itfac && (*itfac)->getMappingType() != Idio && (*itfac)->getMappingType() != NotMapped && (*itfac)->getMappingType() != Proxied) 
			mapped = (*itfac)->mapIt(params, f_PCA_Eqt_Spot, f_PCA_FX_Spot, f_PCA_Eqt_Vol, f_PCA_FX_Vol, f_PCA_IR_Curve, f_PCA_FX_Forward, 
				dts_common, reg_2_cur_map);
	}
}

void FactSet::output_factors_info(std::string dir_output) {
	
	std::ofstream log_file;
	std::string hist_fact_info_fn("hist_fact_info.csv");
	std::string fn(dir_output + hist_fact_info_fn);
    log_file.open(fn.c_str());
	unsigned int i;

	log_file << "Number, Factor name, Mapping type, Factor type, Currency, Region, Number of dates, Number of terms, vol, Resid vol, R2, Regr coef" << std::endl;
	for(i = 0; i < _pos.size(); ++i) {
		log_file << i << ",";
		_pos[i]->print_factor_info(log_file);
	}
	
	log_file.close();
}
void FactSet::output_mapping_info(std::string dir_output) {
	unsigned int i, j, k;
	std::ofstream map_file;
	std::string fn(dir_output + "mapped_factors_info.csv");
	map_file.open(fn.c_str());
	for(i = 0; i < _pos.size(); ++i) {
		if(_pos[i]->getMappingType() == Regressed) {
			map_file << (*_pos[i]).getHistName() << ",";
			std::vector<boost::shared_ptr<Factor> > f_regressed_against((*_pos[i]).getMappedTo());
			for(j = 0; j < f_regressed_against.size(); ++j) 
				map_file << (*f_regressed_against[j]).getHistName() << ",";
			map_file << std::endl;

			map_file << (*_pos[i]).getHistName() << ", vols of the mapped_to:,";
			std::vector<unsigned int> terms((*_pos[i])._mapped_terms);
			for(j = 0; j < f_regressed_against.size(); ++j) {
				std::vector<double> vol((*f_regressed_against[j]).getVol());
				if(terms[j] < vol.size())
					map_file << vol[terms[j]] << ",";
			}
			map_file << std::endl;

			std::vector<double> vol((*_pos[i]).getVol());
			map_file << (*_pos[i]).getHistName() << ", vols:,";
			for(j = 0; j < vol.size(); ++j) 
				map_file << vol[j] << ",";
			map_file << std::endl;

			map_file << (*_pos[i]).getHistName() << ", R2:," << (*_pos[i]).getR2() << std::endl;
			map_file << (*_pos[i]).getHistName() << ", ResidVols:," << (*_pos[i]).getResidVol() << std::endl;

			matrix<double>regrCoef((*_pos[i]).getRegrCoef());
			map_file << (*_pos[i]).getHistName() << ", Regr Coef:,";
			for(j = 0; j < regrCoef.size1(); ++j)  {
				for(k = 0; k < regrCoef.size2(); ++k)
					map_file << regrCoef(j,k) << ",";
			}
			map_file << std::endl;
		}
	}
	map_file.close();
}
int FactSet::output_scenarios2(std::string dir, std::string asOfDate, unsigned int nsim, bool output_diff) {

	std::vector<unsigned int> terms;
	std::vector<double> eod_val;
	matrix<double> rets;
	std::ofstream scen_file;
	double prob(1.0/nsim);
	std::string s_name;
	unsigned int i, j, k, nterms;
	std::ostringstream ostr;

	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::string scen_color("green"), scenSetName("msmc 2000"), scen_name("Mmc2000"), TimeEvolution("Constant"), TimeEvolutionToTrigger("Constant"),
		TriggerHolder("@Standard()"), ScenShiftRule("Trigger Time"), ScenType("non-parallel shift"), ScenReplVal("Term");

	std::vector<boost::shared_ptr<Factor> >::iterator it, jt;
	std::string fn(dir + "scen_MC.csv");
    scen_file.open(fn.c_str());
	
	scen_file << "Scenario Set,Scenario Set Name,Scenario Name,Scenario Probability,Scenario Color,Scenario Variable,Scenario Start Time,Scenario Attribute,Time Evolution,Time Evolution to Trigger,Trigger Holder,Scenario Shift Rule,Scenario Type,Scenario Replacement Value" << std::endl;

	for(i = 0; i < nsim; ++i) {
		std::stringstream out_file;

		j = 0;
		std::ostringstream ostr;
		ostr << (i + 1);
		for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {

			if(itfac->get() == NULL || (*itfac)->getExcludeFromScenarios() == true)
				continue;

			if(i == 0 && j == 0) 
				out_file << ",msmc 2000,";
			else 
				out_file << ",,";

			if(j == 0)
				out_file << "Mmc2000" + ostr.str() << "," << prob << ",";
			else
				out_file << ",,";
			terms = (*itfac)->getTerms();
			nterms = terms.size();
			s_name = (*itfac)->getScenName();

			rets = (*itfac)->_sim_rets;

			out_file << scen_color << "," << s_name << "," << asOfDate << ",,Constant,Constant,@Standard(), Trigger Time,non-parallel shift,Term," << asOfDate << std::endl;
			if(output_diff) {
				eod_val = (*itfac)->getEODVal();
				if(eod_val.size() != nterms)
					for(k = 0; k < nterms; ++k) 
						out_file << ",,,,,,,,,,,,," << terms[k] << "," << 0.0 << std::endl;
				else
					for(k = 0; k < nterms; ++k) 
						if((*itfac)->getType() == FX_Pair_Vol )
							out_file << ",,,,,,,,,,,,," << terms[k] << "," << exp(rets(i,k)) - 1.0 << std::endl;
						else {
							if((*itfac)->isLogNormal())
								out_file << ",,,,,,,,,,,,," << terms[k] << "," << eod_val[k] * (exp(rets(i,k)) - 1.0) << std::endl;
							else
								out_file << ",,,,,,,,,,,,," << terms[k] << "," << rets(i,k) << std::endl;
						}
			}
			else
				for(k = 0; k < nterms; ++k) 
					out_file << ",,,,,,,,,,,,," << terms[k] << "," << rets(i, k) << std::endl;
			++j;
		}
		// scenValue, scenNum
		out_file << ",,,," << scen_color << ",:scenValue," << asOfDate << ",,Constant,Constant,@Standard(), Trigger Time,non-parallel shift,Term," << asOfDate << std::endl;
		out_file << ",,,,,,,,,,,,,0.0,3" << std::endl;

		out_file << ",,,," << scen_color << ",:scenNum," << asOfDate << ",,Constant,Constant,@Standard(), Trigger Time,non-parallel shift,Term," << asOfDate << std::endl;
		out_file << ",,,,,,,,,,,,,0.0," << i+1 << std::endl;

		scen_file << out_file.rdbuf();
	}
	// Base scenario
	j = 0;
	std::ostringstream ostr2;
	ostr2 << (i + 1);
	std::stringstream out_file;
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		if(itfac->get() == NULL || (*itfac)->getExcludeFromScenarios() == true)
				continue;
		if(j == 0)
			out_file << ",,Mmc2000" + ostr2.str() << "," << prob << ",";
		else
			out_file << ",,,,";
		terms = (*itfac)->getTerms();
		nterms = terms.size();
		s_name = (*itfac)->getScenName();

		out_file << scen_color << "," << s_name << "," << asOfDate << ",,Constant,Constant,@Standard(), Trigger Time,non-parallel shift,Term," << asOfDate << std::endl;
			
		for(k = 0; k < nterms; ++k) 
			out_file << ",,,,,,,,,,,,," << terms[k] << "," << 0.0 << std::endl;
		++j;
	}
	// scenValue, scenNum
	out_file << ",,,," << scen_color << ",:scenValue," << asOfDate << ",,Constant,Constant,@Standard(), Trigger Time,non-parallel shift,Term," << asOfDate << std::endl;
	out_file << ",,,,,,,,,,,,,0.0,0" << std::endl;

	out_file << ",,,," << scen_color << ",:scenNum," << asOfDate << ",,Constant,Constant,@Standard(), Trigger Time,non-parallel shift,Term," << asOfDate << std::endl;
	out_file << ",,,,,,,,,,,,,0.0,0" << std::endl;

	scen_file << out_file.rdbuf();
	scen_file.close();
 	return 0;
}

void FactSet::assignEmpirDistFlag(std::map<std::string, std::vector<unsigned int> >& empirDistMap) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::map<std::string, std::vector<unsigned int> >::iterator itmap;
	std::string h_name;
	unsigned int i;
	std::vector<unsigned int> terms;
	std::vector<unsigned int>::iterator itvec;
	std::vector<bool> dist_flags;
	std::vector<unsigned int> terms_dist;

	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		h_name = (*itfac)->getHistName();
		itmap = empirDistMap.find(h_name);
		if(itmap != empirDistMap.end()) {
			terms_dist = itmap->second;
			terms = (*itfac)->getTerms();

			dist_flags = (*itfac)->getUseEmpirDist();
			if(dist_flags.size() != terms.size()) 
				dist_flags.resize(terms.size());
			
			for(i = 0; i < terms_dist.size();  ++i) {
				itvec = find(terms.begin(), terms.end(), terms_dist[i]);
				if(itvec != terms.end()) 
					dist_flags[itvec - terms.begin()] = true;
			}
			(*itfac)->setDistFlags(dist_flags);
		}
	}
}
void FactSet::getFactorsByMappingType(FactMappingType type, std::vector<boost::shared_ptr<Factor> >& o_facs) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		if((*itfac)->getMappingType() == type )
			o_facs.push_back(*itfac);
	}
}
void FactSet::getFactorsByType(FactType factType, std::vector<boost::shared_ptr<Factor> >& f_vec) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		if((*itfac)->getType() == factType && (*itfac)->getNterms() > 0)
			f_vec.push_back(*itfac);
	}		
}
void FactSet::addPositions(std::vector<boost::shared_ptr<Factor> >& facs) {
	for(unsigned int i = 0; i < facs.size(); ++i)
		_pos.push_back(facs[i]);
}
void FactSet::output_hr_factors(std::string fn, unsigned int sim_horizon, double decay) {
	std::ofstream ff;
	std::vector<date> dts(_pos[0]->getHrDates());
	std::string name("");
	unsigned int i(0), j(0);
	double sim_horizon_days(1.0);
	std::vector<double> vals(dts.size()), wts;
	date_duration dt1, dt2;
	matrix<double> hr;
	std::string yyyymmdd, ret;
	char date[20];
	ff.open(fn.c_str());
	calculate_weights(decay, dts, wts, sim_horizon);

	for(i = 0; i < _pos.size(); ++i) {
		ff << "FactorGroup,Name,Currency,Price/Yield,Date,Value" << std::endl;
		std::ostringstream result;

		result << i;
		if( i < 10)
			name = "0"  + result.str();
		else
			name = result.str();
		
		ff << "equity_data,J" << name << "IRVol,J" << name << ",Price,";
		hr = _pos[i]->_rets;
		vals[0] = 1.0;
		sprintf ( date, "%d/%d/%d", static_cast<short>(dts[0].year()), static_cast<short>(dts[0].month()), static_cast<short>(dts[0].day()));
		ff << date << "," << vals[0] << std::endl;
		
		for(j = 1; j < dts.size(); ++j) {
			dt1 = dts[sim_horizon] - dts[0];
			dt2 = dts[j] - dts[0];
			sim_horizon_days = dt1.days();
			if(j < sim_horizon)
				vals[j] = exp(hr(0, 0) * wts[j] * dt2.days() / sim_horizon_days);
			else
				vals[j] = vals[j - sim_horizon] * wts[j - sim_horizon] * exp(hr(j - sim_horizon, 0));
			sprintf (date, "%d/%d/%d", static_cast<short>(dts[j].year()), static_cast<short>(dts[j].month()), static_cast<short>(dts[j].day()));
			ff << ",,,," << date << "," << vals[j] << std::endl;
		}
	}

	ff.close();
}

void FactSet::readProxyMap(std::string fn, double defaultIdioVol, double defaultIdioEquityVol) {
	// This is implemented here for testing only (to compare with current production). These mappings do not make any sense Tanya
	// The proxy name and multiplier are assumed the same for all tenors
	std::ifstream ff;
	
	// That's how the file looks like 
	// For unexplicable reason FX forward curves are mapped to FX spot and FX spot names are represented as in historical data file (ARSSpot), not as all others here (ALGO format)
	// Risk Factor Type, Risk Factor Name, Term, Proxy Name, Proxy Multiplier, Proxy Spread, Comments
	// IR, FXARS_Forward,,ARSSpot, ,,,
	// IR, FXBBD_Forward,,BBDSpot, ,,,

	std::string buf(""), trimmed_name;
	std::vector<std::string> values;
	std::vector<double> scalers, scalers_cur, interpolated_scalers, one(1), vol, eod, arg;
	std::vector<unsigned int> tenors, terms, jterms, proxy_tenors;
	std::string prev_name(""), name, s_name, proxy_name;
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac, jtfac;
	unsigned int nterms(1), i;
	double val(1.0);
	std::map<double, double> mp;

	one[0] = 1.0;
	using boost::is_any_of;

	ff.open(fn.c_str());
	if(ff.is_open()) {
	getline(ff, buf); // First line is headers
	while(getline(ff, buf)) {
		boost::algorithm::split(values, buf, is_any_of(",")); 
		name = values[1];
		trim(name);
		if(name == "IRUSD_Supranational")
			int tt = 0;

		s_name = ":" + name;
		itfac = matchNameWithFactor(s_name, _pos, false);
		if(itfac != _pos.end()) {
			if((*itfac)->getHistName() == "AUDBBS3M" || (*itfac)->getHistName() == "C78462F10_NYSE_SPY_ETFSpot")
				int tt = 0;

			// This values[3] is the proxy name. It is h_name for names that corresond to FX_Spot, like "NADSpot", but it is s_name without ":" for name that correspond to IR like IRAUD_Interbank
			// Make proxy name in all cases to be s_name
			(*itfac)->setMappingType(Proxied);

			trimmed_name = values[3];
			trim(trimmed_name);
			if(trimmed_name.substr(trimmed_name.length() - 4, trimmed_name.length()) == "Spot")
				proxy_name = ":FX" + trimmed_name.substr(0, 3);
			else
				proxy_name = ":" + trimmed_name;
			(*itfac)->setProxyName(proxy_name);
			if(values.size() == 4) { // Example: IR, IRNZD_Interbank, , NZDForward
				(*itfac)->getProxyTenors().push_back(0);
				(*itfac)->getProxyScalers().push_back(1.0);
			}
			if(values.size() >= 7 ) { // Example: IR, IRGBP_RPI, 30, GBPInterbank, -0.239148
				(*itfac)->getProxyTenors().push_back(atoi(values[2].c_str()));
				if(values[4] == "" || values[4] == " ")
					val = 1.0;
				else
					val = atof(values[4].c_str());
				(*itfac)->getProxyScalers().push_back(val);
				if(values[6] == "UseEOD") 
					(*itfac)->setUseProxyEOD(true);
			}
		}
	}
	ff.close();
	}
	// What does proxiing affect
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		jtfac = findProxyFact((*itfac)->getProxyName());
		nterms = (*itfac)->getNterms();
		
		if((*itfac)->getHistName() == "C08467070_NYSE_BRKBSpot" || (*itfac)->getHistName() == "AUDBBS3M")
			int tt = 0;

		if((nterms == 0 && jtfac == _pos.end()) || (nterms == 0 && jtfac != _pos.end() && ((*itfac)->getType() != (*jtfac)->getType() || (*jtfac)->getTerms().size() == 0))) { 
			// No historical data, no proxy This factor will have zero risk because there is no EOD data for now Tanya
			(*itfac)->handle_noData_noProxy(defaultIdioVol, defaultIdioEquityVol);
			continue;
		}
		if(jtfac == _pos.end())
			continue;

		// New (proxy) tenors
		jterms = (*jtfac)->getTerms();
		if(jterms.size() == 0)
			continue;

		arg.resize(jterms.size());
		for(i = 0; i < jterms.size(); ++i) 
			arg[i] = jterms[i];

		// Set EOD from historical data for the new (proxy) tenors. EOD will change later if there is EOD data
		mp.clear();
		if(nterms > 0) {// Set this factors' EOD values at proxy tenors
			for(i = 0; i < nterms; ++i) 
				mp.insert(std::make_pair<double,double>((*itfac)->getTerms()[i], (*itfac)->getEODVal()[i]));
			
			interpolate_linear(mp, arg, eod);	
		}
		else  // Set proxy's EOD values at proxy tenors
			eod = (*jtfac)->getEODVal();

		(*itfac)->setEODVal(eod);

		// Set vols, terms and proxy scalers as in The Proxy
		mp.clear();
		scalers_cur = (*itfac)->getProxyScalers();
		proxy_tenors = (*itfac)->getProxyTenors();
		for(i = 0; i < scalers_cur.size(); ++i) 
			mp.insert(std::make_pair<double, double>(proxy_tenors[i], scalers_cur[i]));

		interpolate_linear(mp, arg, scalers);
		(*itfac)->setProxyScalers(scalers);
			
		(*itfac)->setTerms((*jtfac)->getTerms());
		(*itfac)->setVol((*jtfac)->getVol());
			
			/* 
			std::vector<double> vols, proxy_terms_dbl, factor_terms_dbl;
			std::vector<unsigned int> proxy_terms((*itfac)->getProxyTenors()), factor_terms((*itfac)->getTerms());

			for(i = 0; i < proxy_terms.size(); ++i) 
				proxy_terms_dbl.push_back(proxy_terms[i]);

			for(i = 0; i < factor_terms.size(); ++i)
				factor_terms_dbl.push_back(factor_terms[i]);
			
			interpolate_linear(proxy_terms_dbl, (*itfac)->getProxyScalers(), factor_terms_dbl, interpolated_scalers);
			(*itfac)->setProxyScalers(interpolated_scalers);

			proxy_terms_dbl.clear();
			proxy_terms = (*jtfac)->getTerms();
			for(i = 0; i < proxy_terms.size(); ++i) 
				proxy_terms_dbl.push_back(proxy_terms[i]);

			interpolate_linear(proxy_terms_dbl, (*jtfac)->getVol(), factor_terms_dbl, vols);
			(*itfac)->setVol(vols);
		*/
	}
}
/*
std::vector<boost::shared_ptr<Factor> >::iterator FactSet::matchNameWithFactor(std::string name, bool isHistName) {

	std::vector<boost::shared_ptr<Factor> >::iterator itfac(_pos.begin());
	if(isHistName)
		while(itfac != _pos.end() && name != (*itfac)->getHistName())
			++itfac;
	else
		while(itfac != _pos.end() && name != (*itfac)->getScenName())
			++itfac;
	return itfac;
}
*/
std::vector<boost::shared_ptr<Factor> >::const_iterator FactSet::findProxyFact(std::string proxy_name) const {
	std::vector<boost::shared_ptr<Factor> >::const_iterator jtfac;
	jtfac = matchNameWithFactor(proxy_name, _pos, false);
	return jtfac;
}
void FactSet::genSimsProxied(std::vector<boost::shared_ptr<Factor> >& proxied, unsigned int nsim) {
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac, jtfac;
	std::string name("");
	std::vector<unsigned int> terms, jterms;
	std::vector<unsigned int>::iterator it1, it2;
	std::vector<double> scalers, scalers_cur, arg, eod, eod_cur;
	std::vector<unsigned int> proxy_tenors;
	std::map<double, double> mp;

	unsigned int nterms, isim, i;
	double w1(0.0), w2(0.0), x1(0.0), x2(0.0), x(0.0);

	for(itfac = proxied.begin(); itfac != proxied.end(); ++itfac) {
		nterms = (*itfac)->getNterms();
		(*itfac)->_sim_rets.resize(nsim, nterms);
		(*itfac)->_sim_rets.clear();
		jtfac = findProxyFact((*itfac)->getProxyName());

		if((*itfac)->getHistName() == "EUREIB1Y")
			int tt = 0;

		if(jtfac == _pos.end()) 
			continue; // Proxy factor is not found

		terms = (*itfac)->getTerms();
		jterms =(*jtfac)->getTerms();

		if(jterms.size() == 0) {
			for(i = 0; i < nterms; ++i)
				for(isim = 0; isim < nsim; ++isim)
					(*itfac)->_sim_rets(isim, i) = 0.0;
				
				continue;
		}
		// In current production tenors of the proxy are used
		scalers = (*itfac)->getProxyScalers();
		(*itfac)->_sim_rets.resize((*jtfac)->_sim_rets.size1(), (*jtfac)->_sim_rets.size2());
		for(i = 0; i < jterms.size(); ++i)
			column((*itfac)->_sim_rets, i) = column((*jtfac)->_sim_rets, i) * scalers[i];
		/*
		eod.resize(jterms.size());
		eod_cur = (*itfac)->getEODVal();	
		mp.clear();
		for(i = 0; i < nterms; ++i) 
			mp.insert(std::make_pair<double, double>(terms[i], eod_cur[i]));
			
		interpolate_linear(mp, arg, eod);
		(*itfac)->setTerms((*jtfac)->getTerms());
		(*itfac)->setVol((*jtfac)->getVol());
		(*itfac)->setEODVal(eod);
		*/
		/* In this version the tenors for the proxied factor were used
		if(jterms.size() == 1) {
			for(i = 0; i < nterms; ++i)
				column((*itfac)->_sim_rets, i) = column((*jtfac)->_sim_rets, 0);

			continue;
		}
		for(i = 0; i < nterms; ++i) {
				it1 = lower_bound(jterms.begin(), jterms.end(), terms[i]);
				it2 = upper_bound(jterms.begin(), jterms.end(), terms[i]);
	
				if(it2 == jterms.begin() ) 
					column((*itfac)->_sim_rets, i) = column((*jtfac)->_sim_rets, 0); 
				else if (it1 == jterms.end())
					column((*itfac)->_sim_rets, i) = column((*jtfac)->_sim_rets, jterms.size() - 1);
				else {
					if(it2 == jterms.end())
						i2 = jterms.size() - 1;
					else
						i2 = it2 - jterms.begin();
					i1 = i2 - 1;	
				}
				x2 = jterms[i2];
				x1 = jterms[i1];
				x = terms[i];
				w1 = (x2 - x) / (x2 - x1);
				w2 = 1.0 - w1;
				column((*itfac)->_sim_rets, i) = column((*jtfac)->_sim_rets, i1) * w1 + column((*jtfac)->_sim_rets, i2) * w2;
		}
		*/
	}
}
void FactSet::removeEqtFactors() {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;

	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		if((*itfac)->getType() == Eqt_Spot && (*itfac)->getMappingType() != PCA) 
			(*itfac)->setExcludeFromScenarios(true);
	}
}
void FactSet::calc_factor_coef(std::map<unsigned int, std::string>& stockMap, std::map<unsigned int, std::vector<double> >& factorCoefMap, 
		std::vector<boost::shared_ptr<Factor> >& pca, unsigned int sim_horizon, double decay) {

	unsigned num(0), j(0), nrets(0), i, npca(pca.size());
	matrix<double> regr_coef;
	std::vector<double> vec(npca + 1), wts;
	std::map<std::string, std::vector<double> > regrCoefMap;
	double sum(0.0), vol_resid(0.0), systematic_var(0.0), pca_vol(0.0), xTy(0.0), vol(0.0);
	std::string name;
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::map<unsigned int, std::string>::iterator itStockMap;
	std::map<std::string, std::vector<double> >::iterator itRegrCoef;
	std::map<unsigned int, std::vector<double> >::iterator itFactorCoef;
	boost::numeric::ublas::vector<double> Y, X;

	if(npca == 0) {
		for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
			name = (*itfac)->getHistName();
			vec[0] = (*itfac)->getVol()[0];
			regrCoefMap.insert(std::make_pair(name, vec));
		}
	}
	else {
		factorCoefMap.clear();
		calculate_weights(decay, pca[0]->getHrDates(), wts, sim_horizon);
		nrets = wts.size();
		Y.resize(nrets);
		X. resize(nrets);
		regr_coef.resize(npca, 1);

		for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {

			eig_info << (*itfac)->getHistName() << ",";
			if((*itfac)->getType() == Eqt_Spot) {
				if((*itfac)->getMappingType() == NotMapped) {
					for(i = 0; i < nrets; ++i)
						Y[i] = (*itfac)->_rets(i, 0) * wts[i];
					vol = (*itfac)->getVol()[0];
					systematic_var = 0.0;
					for(j = 0; j < npca; ++j) {
						for(i = 0; i < nrets; ++i)
							X[i] = pca[j]->_rets(i, 0) * wts[i];
						pca_vol = pca[j]->getVol()[0];

						xTy = inner_product(Y.begin(), Y.end(), X.begin(), 0.0);
						regr_coef(j, 0) = xTy / pca_vol;
						systematic_var += xTy * xTy;
					}
					
					vol_resid = vol * vol - systematic_var;
					if(vol_resid > 0)
						vol_resid = sqrt(vol_resid);

					(*itfac)->setRegrCoef(regr_coef);
					(*itfac)->setResidVol(vol_resid);
				}
				if((*itfac)->getMappingType() == Regressed || (*itfac)->getMappingType() == NotMapped) {
					regr_coef = (*itfac)->getRegrCoef();
					if(regr_coef.size1() != npca) 
						regr_coef.resize(npca, 0.0);

					vec[0] = (*itfac)->getResidVol();
					for(j = 1; j < vec.size(); ++j)
						vec[j] = regr_coef(j - 1, 0);
				}
				else if((*itfac)->getMappingType() == Idio) {
					for(j = 1; j < vec.size(); ++j)
						vec[j] = 0.0;
					vec[0] = (*itfac)->getVol()[0];
				}
				else
					continue;
			
				name = (*itfac)->getHistName();
				regrCoefMap.insert(std::make_pair(name, vec));
				eig_info << name << "," << vec.size() << std::endl;
			}
		}
	}
	/*
	eig_info << "factorCoefMap.insert" << std::endl;
	for(itRegrCoef = regrCoefMap.begin(); itRegrCoef != regrCoefMap.end(); ++itRegrCoef) {
		name = itRegrCoef->first;
		itStockMap = stockMap
	}
	*/
	for(itStockMap = stockMap.begin(); itStockMap != stockMap.end(); ++itStockMap) {
		num = itStockMap->first;
		name = itStockMap->second;
		eig_info << num << "," << name << ",";
		itRegrCoef = regrCoefMap.find(name);
		if(itRegrCoef != regrCoefMap.end()) {
			vec = itRegrCoef->second;
			factorCoefMap.insert(std::make_pair(num, vec));
			eig_info << vec.size() << "got it,";
		}
		eig_info << std::endl;
	}
}
void FactSet::output_factor_coef(std::map<unsigned int, std::string>& stockMap, std::map<unsigned int, std::vector<double> >& factorCoefMap, std::string fn) {
	
	std::ofstream ff;
	std::map<unsigned int, std::string>::iterator itStockMap;
	std::map<unsigned int, std::vector<double> >::iterator itFactorCoef;
	std::string name;
	std::vector<double> vec;
	unsigned int num(0), j;

	ff.open(fn.c_str());
	ff << "Matrix Parameter Curve, Name, ID, Datum, RelativeCurve,Time Evolution,Interpolate Axis 1, Interpolate Axis 2, CurveUnit,Surface,,,,," << std::endl;
	ff << ",IndexCoef_GB,IndexCoef_GB,0,True,Constant,Constant,Constant,% ANNU actual/actual,Axis 1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18," << std::endl;
	
	for(itFactorCoef = factorCoefMap.begin(); itFactorCoef != factorCoefMap.end(); ++itFactorCoef) {
		num = itFactorCoef->first;
		itStockMap = stockMap.find(num);
		name = itStockMap->second;
		//ff << ",,,,,,,,," << itFactorCoef->first;
		ff << ",,,,,,,," << name << "," << itFactorCoef->first;
		vec = itFactorCoef->second;
		for(j = 0; j < vec.size(); ++j)
			ff << "," << vec[j];
		ff << std::endl;
	}
	ff.close();
}
void FactSet::setEOD_ProxiedUseProxyEOD() {
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac, jtfac;
	unsigned int iterm;
	std::map<double, double> mp;
	std::vector<double> val, arg, eod_val;

	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		if(!(*itfac)->getUseProxyEOD())
			continue;
		
		jtfac = findProxyFact((*itfac)->getProxyName());
		if(jtfac != getFactors().end() && (*itfac)->getType() == (*jtfac)->getType()) {
			arg.resize((*jtfac)->getNterms());
			val = (*jtfac)->getEODVal();
			for(iterm = 0; iterm < arg.size(); ++iterm)
				mp.insert(std::make_pair<double, double>(arg[iterm], val[iterm]));

			interpolate_fwd_linear(mp, arg, eod_val);
			(*itfac)->setEODVal();
		}
	}
}
void FactSet::genSensShifts_all(std::string asOfDate, std::string dir, double shift, unsigned int npca) {
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac;
	std::string ir_fn("shifts_ir.csv"), non_ir_fn("shifts_non_ir.csv"), ir_info_fn("shifts_ir_info.csv"), non_ir_info_fn("shifts_non_ir_info.csv");
	bool doIR(true);
	genShifts(dir + ir_fn, dir + ir_info_fn, asOfDate, shift, "SensDelta", doIR, true, true);

	// When writing shifts we have to comply with current naming convention for :factorGB_0, etc, which are assumed to be 18 different factors
	// When writing shift info file we switch to have a single :factorGB risk factor with 18 tenors. Then process to calculate gammas is the same at for all other factors with multiple tenors
	genShifts(dir + non_ir_fn, dir + non_ir_info_fn, asOfDate, shift, "SensDelta", !doIR, true, false);

	// Put PCA in one factor with multiple tenors
	makeSinglePCAFactor(npca);
	genShifts(dir + non_ir_fn, dir + non_ir_info_fn, asOfDate, shift, "SensDelta", !doIR, false, true);
}

void FactSet::apply_CCAR_info(const std::vector<boost::shared_ptr<CCAR_shift> >& ccar_shifts, unsigned int i_proj, 
	const std::map<std::pair<std::string, std::string>, std::string>& mp, const std::map<std::string, std::string>& reg_2_cur_mp, const std::map<std::string, std::string>& cur_2_reg_map) {
							
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac;
	
	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) 
		(*itfac)->apply_CCAR_info(ccar_shifts, i_proj, mp, reg_2_cur_mp, cur_2_reg_map);
}
int FactSet::genShifts(std::string file_nm, std::string factors_info_fn, std::string asOfDate, double shift, std::string scenName, bool doIR, bool genShifts, bool getShiftsInfo) {
	// Shifts will be written to a file plus info in  fact_info file
	std::vector<unsigned int> iterms, jterms;
	std::ofstream fact_info_file, scen_file;
	
	const std::string line1("Scenario Set,Scenario Set Name,Scenario Name,Scenario Probability,Scenario Color,Scenario Variable,Scenario Start Time,Scenario Attribute,Time Evolution,Time Evolution to Trigger,Trigger Holder,	Scenario Shift Rule,Scenario Type,Scenario Replacement Value,");
	unsigned int i, j, k(1), count;
	double prob(0.001), this_shift;
	bool checkType;
	const std::string scen_color("green"), scenSetName("msmc 2000"), scen_name(scenName + "000"), ScenarioStartTime(asOfDate), TimeEvolution("Constant"), TimeEvolutionToTrigger("Constant"),
		TriggerHolder("@Standard()"), ScenShiftRule("Trigger Time"), ScenType("non-parallel shift"), ScenReplVal("Term"), Attr14(asOfDate);
	std::string factName1(""), factName2("");

	std::vector<boost::shared_ptr<Factor> >::iterator it, jt;
	FactType iType;
	ShiftType shiftType(DeltaUp);
	iType = _pos[0]->getType();

	if(genShifts) {
		scen_file.open(file_nm);
		scen_file << line1 << std::endl;
	}
	if(getShiftsInfo) {
		fact_info_file.open(factors_info_fn);
		fact_info_file << "ShiftType, Scenario, Factor 1, Factor 2, tenor1, tenor2, shift" << std::endl;
	}

	// Delta scenarios
	k = 1;
	for(count = 0; count < 2; ++count) {
		switch(count) {
		case 0:
			shiftType = DeltaDn;
			break;
		case 1:
			shiftType = DeltaUp;
			break;
		}
		for(it = _pos.begin(); it != _pos.end(); ++it) {
			
			checkType = (*it)->getType() == IR_Curve || (*it)->getType() == IRCap_Vol;
			if(doIR && !checkType || !doIR && checkType) 
				continue;

			if((*it)->getType() == Eqt_Spot && (*it)->getMappingType() != PCA)
				continue;

			factName1 = (*it)->getScenName();
			iterms = (*it)->getTerms();
			for(i = 0; i < iterms.size(); ++i) { // new scenario
				
				//this_shift = shift * (*it)->getEODVal()[0];
				this_shift = shift;
				if(genShifts) {
					std::stringstream out_file;
					output_scenario(scenSetName, scen_name, scen_color, ScenarioStartTime, TimeEvolution, TimeEvolutionToTrigger, TriggerHolder, ScenShiftRule, ScenType,ScenReplVal, Attr14, 
						shiftType, k, prob, this_shift, doIR, it, it, i, i, out_file);
					scen_file << out_file.rdbuf();
				}
				if(getShiftsInfo)
					fact_info_file << shiftType << "," << k-1 << "," << factName1 << "," << factName1 << "," << iterms[i] << "," << iterms[i] << "," << this_shift << std::endl;
				
				++k;
			} // end for
		}
	}
	// Gamma scenarios
	for(count = 0; count < 4; ++count) {
		switch(count) {
		case 0:
			shiftType = CrossGammaUpUp;
			break;
		case 1:
			shiftType = CrossGammaUpDn;
			break;
		case 2:
			shiftType = CrossGammaDnUp;
			break;
		case 3:
			shiftType = CrossGammaDnDn;
			break;
		}
		for(it = _pos.begin(); it != _pos.end(); ++it) {
			checkType = (*it)->getType() == IR_Curve || (*it)->getType() == IRCap_Vol;
			if(doIR && !checkType || !doIR && checkType) 
				continue;

			if((*it)->getType() == Eqt_Spot && (*it)->getMappingType() != PCA)
				continue;

			factName1 = (*it)->getScenName();
			factName2 = factName1;
			
			iterms = (*it)->getTerms();
			for(i = 0; i < iterms.size(); ++i) {
				
				for(j = i + 1; j < iterms.size(); ++j) {  // new scenario
					
					//this_shift = shift * (*it)->getEODVal()[0];
					this_shift = shift;
					if(genShifts) {
						std::stringstream out_file;
						output_scenario(scenSetName, scen_name, scen_color, ScenarioStartTime, TimeEvolution, TimeEvolutionToTrigger, TriggerHolder, ScenShiftRule, ScenType,ScenReplVal, Attr14, 
							shiftType, k, prob, this_shift, doIR, it, it, i, j, out_file);
						scen_file << out_file.rdbuf();
					}
					if(getShiftsInfo) 
						fact_info_file << shiftType << "," << k-1 << "," << factName1 << "," << factName2 << "," << iterms[i] << "," << iterms[j] << "," << this_shift << std::endl;
					
					++k;
				}
			}
		}
		if(!doIR && genShifts) {  // This is for Equity PCA now
			i = 0;
			j = 0;
			for(it = _pos.begin(); it != _pos.end(); ++it) {
				if((*it)->getMappingType() == PCA) {
					for(jt = it + 1; jt != _pos.end(); ++jt) {
						if((*jt)->getMappingType() == PCA) { // new scenario
							
								std::stringstream out_file;
								output_scenario(scenSetName, scen_name, scen_color, ScenarioStartTime, TimeEvolution, TimeEvolutionToTrigger, TriggerHolder, ScenShiftRule, ScenType,ScenReplVal, Attr14, 
									shiftType, k, prob, shift, doIR, it, jt, i, j, out_file);
								scen_file << out_file.rdbuf();
							
							//	fact_info_file << shiftType << "," << k-1 << "," << ":factorGB,:factorGB," << i << "," << j << shift << std::endl;
				
							++j;
							++k;
						}
					}
					++i;
				}
			}
		}
	}			
	// Add base scenario
	if(genShifts) {
		std::stringstream out_file;
		output_scenario(scenSetName, "Base", scen_color, ScenarioStartTime, TimeEvolution, TimeEvolutionToTrigger, TriggerHolder, ScenShiftRule, ScenType,ScenReplVal, Attr14, 
			Base, 0, prob, shift, doIR, it, it, 0, 0, out_file);
		scen_file << out_file.rdbuf();
		scen_file.close();
	}
	if(getShiftsInfo) {
		fact_info_file << "0," << k-1 << ",Base, Base, 0, 0, 0.0" << std::endl;
		fact_info_file.close();
	}
 	return 0;
}

void FactSet::output_scenario(std::string scenSetName, std::string scen_name, std::string scen_color, std::string ScenarioStartTime, std::string TimeEvolution, std::string TimeEvolutionToTrigger,
	std::string TriggerHolder, std::string ScenShiftRule, std::string  ScenType,std::string  ScenReplVal, std::string asOfDate, 
		ShiftType shiftType, unsigned int k, double prob, double shift, bool doIR, std::vector<boost::shared_ptr<Factor> >::iterator it, 
		std::vector<boost::shared_ptr<Factor> >::iterator mt, unsigned int i0,  unsigned int j0, std::stringstream & out_file) {
	bool checkType;
	double shiftVal;
	std::vector<boost::shared_ptr<Factor> >::iterator jt;
	unsigned int m(0), j, scenNum(k), scenVal(2);
	std::vector<unsigned int> jterms;
	std::string factName;

	if(k == 1)
		out_file << "," << scenSetName << ",";
	else
		out_file << ",,";

	if(shiftType == Base)
		out_file << "Base," << prob;
	else
		out_file << scen_name << k << "," << prob;
	
	for(jt = _pos.begin(); jt != _pos.end(); ++jt) {

		checkType = (*jt)->getType() == IR_Curve || (*jt)->getType() == IRCap_Vol;
		if(doIR && !checkType || !doIR && checkType) 
				continue;
		if((*jt)->getType() == Eqt_Spot && (*jt)->getMappingType() != PCA)
				continue;

		factName = (*jt)->getScenName();
		out_file << ",";
		if(m != 0)
			out_file << ",,,";
		out_file << scen_color << "," << factName << "," << ScenarioStartTime << ",," << TimeEvolution << "," << TimeEvolutionToTrigger << "," << 
			TriggerHolder << "," << ScenShiftRule << "," << ScenType << "," << ScenReplVal << "," << asOfDate << std::endl;
		jterms = (*jt)->getTerms();
		for(j = 0; j < jterms.size(); ++j) {
			shiftVal = 0.0;
			if(k != 0 && jt == it && i0 == j0 && j == i0) {
				if(shiftType == DeltaUp)
					shiftVal = shift;
				else
					shiftVal = -shift;
			}
			if(k != 0 && it == mt && jt == it && i0 != j0) {
				if(shiftType == CrossGammaUpUp && (j == i0 || j == j0))
					shiftVal = shift;
				if(shiftType == CrossGammaDnDn && (j == i0 || j == j0))
					shiftVal = -shift;
				if((shiftType == CrossGammaUpDn && j == i0) || (shiftType == CrossGammaDnUp && j == j0))
					shiftVal = shift;
				if((shiftType == CrossGammaUpDn && j == j0) || (shiftType == CrossGammaDnUp && j == i0))
					shiftVal = -shift;				
			}
			if(k != 0 && it != mt && (jt == it || jt == mt)) { // For factors with single term
				if(shiftType == CrossGammaUpUp && (jt == it || jt == mt))
					shiftVal = shift;
				if(shiftType == CrossGammaDnDn && (jt == it || jt == mt))
					shiftVal = -shift;
				if((shiftType == CrossGammaUpDn && jt == it) || (shiftType == CrossGammaDnUp && jt == mt))
					shiftVal = shift;
				if((shiftType == CrossGammaUpDn && jt == mt) || (shiftType == CrossGammaDnUp && jt == it))
					shiftVal = -shift;				
			}
			out_file << ",,,,,,,,,,,,," << jterms[j] << "," << shiftVal << std::endl;
		}
		++m;
	} // end for
	// Add scenValue, scenNum
	if(shiftType == Base) {
		scenNum = 0;
		scenVal = 0;
	}
	out_file << ",,,," << scen_color << ",:scenValue," << ScenarioStartTime << ",," << TimeEvolution << "," << TimeEvolutionToTrigger << "," << 
		TriggerHolder << "," << ScenShiftRule << "," << ScenType << "," << ScenReplVal << "," << asOfDate << std::endl;
	out_file << ",,,,,,,,,,,,,0.0," << scenVal << std::endl;

	out_file << ",,,," << scen_color << ",:scenNum," << ScenarioStartTime << ",," << TimeEvolution << "," << TimeEvolutionToTrigger << "," << 
		TriggerHolder << "," << ScenShiftRule << "," << ScenType << "," << ScenReplVal << "," << asOfDate << std::endl;
	out_file << ",,,,,,,,,,,,,0.0," << scenNum  << std::endl;
}

void FactSet::resetSens() {
	unsigned int nterms;
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac;

	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		(*itfac)->_delta.clear();
		(*itfac)->_gamma.clear();
		nterms = (*itfac)->getTerms().size();
		(*itfac)->_delta.resize(nterms, 0.0);
		(*itfac)->_gamma.resize(nterms, nterms);
		(*itfac)->_gamma *= 0.0;
	}
}
bool FactSet::setSens_fromShiftVals(const std::map<unsigned int, ShiftScenInfo>& shift_info_mp, std::map<unsigned int, double>& val_mp, double shift) {
	// Uses shift scenario info to interpret scenarios and calculate sensitivities
	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac;
	std::map<unsigned int, ShiftScenInfo>::const_iterator it_info_mp;
	std::map<unsigned int, double>::const_iterator it_val_mp;
	unsigned int tenor1, tenor2, nterms, i, j;
	std::string name;
	ShiftScenInfo shiftInfo;
	ShiftType shiftType;
	double baseVal, diff, d1, d2, prec(0.0001);
	std::vector<unsigned int> terms;
	std::vector<unsigned int>::iterator it;

	d1 = prec / shift;
	d2 = prec / (shift * shift);

	// baseVal is assumed to be the last scenario
	baseVal = val_mp[val_mp.size() - 1];

	// Loop over val file lines
	for(it_val_mp = val_mp.begin(); it_val_mp != val_mp.end(); ++it_val_mp) {
		it_info_mp = shift_info_mp.find(it_val_mp->first);
		if(it_info_mp == shift_info_mp.end())
			return false;	// this scenario # was not found in info file

		shiftInfo = it_info_mp->second;
		name = shiftInfo.getName1();
		tenor1 = shiftInfo.getTenor1();
		tenor2= shiftInfo.getTenor2();
		shiftType = shiftInfo.getShiftType();
		shift = shiftInfo.getShift();

		if(shiftType == Base)
			continue;

		itfac = findFactor(name, false);
		if(itfac == _pos.end())
			continue;

		terms = (*itfac)->getTerms();
		nterms = terms.size();
		
		i = 0;
		while(i != nterms && terms[i] != tenor1)
			++i;
		if(i == nterms)
			continue;

		diff = it_val_mp->second - baseVal;
		if(it_val_mp->first == 335)
			log_ff << name << "," << it_val_mp->first << "," << baseVal << "," << i << "," << shift << "," << diff << "," << d1 << "," << d2 << "," << (*itfac)->_delta[i] << std::endl;

		if(abs(diff)>0) {
			if(shiftType == DeltaUp) {
				(*itfac)->_delta[i] += 0.5 * diff / shift;
				(*itfac)->_gamma(i, i) += diff / (shift * shift);
			}
			if(shiftType == DeltaDn) {
				(*itfac)->_delta[i] -= 0.5 * diff / shift;
				(*itfac)->_gamma(i, i) += diff / (shift * shift);
			}
			if(fabs((*itfac)->_delta[i]) < d1)
				(*itfac)->_delta[i] = 0.0;
			if(fabs((*itfac)->_gamma(i, i)) < d2) 
				(*itfac)->_gamma(i, i) = 0.0;

			j = 0;
			while(j != nterms && terms[j] != tenor2)
				++j;
			if(j == nterms)
				continue;

			if(shiftType == CrossGammaUpUp || shiftType == CrossGammaDnDn) {
				(*itfac)->_gamma(i, j) += diff / (shift * shift);
				(*itfac)->_gamma(j, i) = (*itfac)->_gamma(i, j);
			}
			if(shiftType == CrossGammaUpDn || shiftType == CrossGammaDnUp) {
				(*itfac)->_gamma(i, j) -= diff / (shift * shift);
				(*itfac)->_gamma(j, i) = (*itfac)->_gamma(i, j);
			}
			if(fabs((*itfac)->_gamma(i, j)) < d2) {
				(*itfac)->_gamma(i, j) = 0.0;
				(*itfac)->_gamma(j, i) = 0.0;
			}
			
		}
	}
	
	return true;
}
void FactSet::makeSinglePCAFactor(unsigned int num_pca) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac(_pos.begin()), jtfac;
	std::string name;
	std::vector<std::string> values;
	std::vector<unsigned int> terms(num_pca);
	std::vector<bool> dist_flags(num_pca, true);
	std::vector<double> eod(num_pca, 1.0), vol(num_pca);
	unsigned int i, num, nsim, nrets, count(0);
	std::vector<boost::shared_ptr<Factor> > new_pos;
	boost::shared_ptr<Factor> pca(new Factor);

	for(itfac = _pos.begin(); itfac != _pos.end(); ++itfac) {
		
		if((*itfac)->getType() == Eqt_Spot && (*itfac)->getMappingType() == PCA) {
			if(count == 0) {
				nrets = (*itfac)->_rets.size1();
				nsim = (*itfac)->_sim_rets.size1();
				pca->setMappingType(PCA);
				pca->setType(Eqt_Spot);
				pca->setNterms(num_pca);
				for(i = 0; i < num_pca; ++i)
					terms[i] = i;
				pca->setTerms(terms);
				pca->setEODVal(eod);
				pca->_rets.resize(nrets, num_pca);
				pca->_sim_rets.resize(nsim, num_pca);
				pca->setScenName(":factorGB");
				pca->setHistName("JIRVol");
				pca->setCurr("USD");
				
			}
			name = (*itfac)->getScenName();
			boost::algorithm::split(values, name, boost::is_any_of("_"));
			if(values.size() < 2)
				continue;
			num = atoi(values[1].c_str());
			column(pca->_rets, num) = column((*itfac)->_rets, 0);
			column(pca->_sim_rets, num) = column((*itfac)->_sim_rets, 0);
			vol[num] = (*itfac)->getVol()[0];
			++count;
		}
		else {
			boost::shared_ptr<Factor> ff(*itfac);
			new_pos.push_back(std::move(ff));
		}
	}
	if(count > 0) {
		pca->setVol(vol);
		pca->setDistFlags(dist_flags);
		new_pos.push_back(std::move(pca));
	}
	_pos.clear();
	addPositions(new_pos);
}
void FactSet::readSims(std::string fn, unsigned int nsim) {
	std::string buf(""), s_name_prev(""), s_name(""), h_name(""), str(""), curr(""), reg("");
	std::vector<std::string> values;
	unsigned int	nterms, i, jterm(0), nn, nHistRets(0);
	bool isLognormal;
	std::ifstream	ff;
	matrix<double>	sims, hist_rets;
	FactMappingType	mapping_type;
	FactType		fact_type;
	std::vector<bool> dist_flags;
	std::vector<unsigned int> terms;
	std::vector<double> eod, vol;
	std::vector<boost::shared_ptr<Factor> >facs;

	ff.open(fn);

	if(!ff.is_open())
		return;

	getline(ff, buf);

	while(getline(ff, buf)) {
		boost::algorithm::split(values, buf, is_any_of(",")); 
		if(values.size() < nsim + 15)
			continue;

		s_name = values[1].c_str();

		if (s_name != s_name_prev) {
			if(s_name_prev != "") { // save previous
				boost::shared_ptr<Factor> ptr(new Factor);

				ptr->setHistName(values[0].c_str());
				ptr->setScenName(s_name_prev);
				ptr->setHistName(h_name);
				ptr->setNterms(nterms);
				ptr->setTerms(terms);
				ptr->setCurr(curr);
				ptr->setRegion(reg);
				ptr->setEODVal(eod);
				ptr->setVol(vol);
				ptr->setDistFlags(dist_flags);
				ptr->setMappingType(mapping_type);
				ptr->setType(fact_type);
				ptr->setIsLognormal(isLognormal);
				ptr->_sim_rets = sims;
				ptr->_rets = hist_rets;
				facs.push_back(std::move(ptr));
			}
			s_name_prev = s_name;
			h_name = values[0].c_str();
			
			str = values[5].c_str();
			if(str == "Idiosyncratic")
				mapping_type = Idio;
			if(str == "Not Mapped")
				mapping_type = NotMapped;
			if(str == "PCA")
				mapping_type = PCA;
			if(str == "Proxied")
				mapping_type = Proxied;
			if(str == "Regressed")
				mapping_type = Regressed;

			str = values[6].c_str();
			if(str == "IR_Curve")
				fact_type = IR_Curve;
			if(str == "EQ_Vol")
				fact_type = EQ_Vol;
			if(str == "Eqt_Spot")
				fact_type = Eqt_Spot;
			if(str == "FX_Forward")
				fact_type = FX_Forward;
			if(str == "FX_Pair_Vol")
				fact_type = FX_Pair_Vol;
			if(str == "FX_Spot")
				fact_type = FX_Spot;
			if(str == "IRCap_Vol")
				fact_type = IRCap_Vol;
			if(str == "MBS")
				fact_type = MBS;

			nn = atoi(values[7].c_str());
			if(nn == 1)
				isLognormal=true;
			else
				isLognormal = false;

			nHistRets = atoi(values[10].c_str());	
			curr = values[3].c_str();
			reg = values[4].c_str();
			nterms = atoi(values[2].c_str());
			terms.resize(nterms);
			sims.resize(nsim, nterms);
			dist_flags.resize(nterms);
			eod.resize(nterms);
			vol.resize(nterms);
			if(nHistRets > 0)
				hist_rets.resize(nHistRets, nterms);
		}
		jterm = atoi(values[8].c_str());
		terms[jterm] = atoi(values[9].c_str());
		
		for(i = 0; i < nsim; ++i)
			sims(i, jterm) =atof(values[i + 15].c_str());
		
		str = values[11].c_str();
		if(str == "Real Dist")
			dist_flags[jterm] = true;
		else
			dist_flags[jterm] = false;
			
		eod[jterm] = atof(values[13].c_str());
		vol[jterm] = atof(values[12].c_str());
	}
	boost::shared_ptr<Factor> ptr1(new Factor);
	ptr1->setHistName(values[0].c_str());
	ptr1->setScenName(s_name);
	ptr1->setHistName(h_name);
	ptr1->setNterms(nterms);
	ptr1->setTerms(terms);
	ptr1->setCurr(curr);
	ptr1->setRegion(reg);
	ptr1->setEODVal(eod);
	ptr1->setVol(vol);
	ptr1->setDistFlags(dist_flags);
	ptr1->setMappingType(mapping_type);
	ptr1->setType(fact_type);
	ptr1->setIsLognormal(isLognormal);
	ptr1->_sim_rets = sims;
	ptr1->_rets = hist_rets;
	facs.push_back(std::move(ptr1));
	addPositions(facs);
}
/*
void FactSet::set_and_eliminate_factType(const FactSet & pfl, FactType factType) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac(_pos.begin()), jtfac;
	while(itfac != pfl._pos.end()) {
		jtfac = itfac + 1;
		if((*itfac)->getType() == factType && (*itfac)->getMappingType() != PCA) {
			_pos.erase(itfac);
			
		}
		itfac = jtfac;
	}
}
*/
/*
std::string FactSet::setDeltasGammasFromShiftVals(std::vector<boost::shared_ptr<Factor> >& factors, std::map<std::string, std::vector<double> >& pfl_vals, double shift) {
	std::map<std::string, std::vector<double> >::iterator itmap(pfl_vals.find(_name));

	if(itmap == pfl_vals.end())
		return false;

	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::vector<double> prc(itmap->second), del, gam;
	int ll((prc.size() - 1) / 2);
	unsigned int iterm, j(0);
	double delta, gamma, base(prc[prc.size() - 1]), scale_delta(1.0 / shift), scale_gamma(0.5 / shift / shift);
	if(ll < 1 || (2 *ll + 1) != prc.size())
		return false;

	for(itfac = factors.begin(); itfac != factors.end(); ++itfac) {
		del.clear();
		gam.clear();
		for(iterm = 0; iterm < (*itfac)->getNterms(); ++iterm) {
			if(j + ll >= prc.size()) 
				return std::string("setDeltasGammasFromShiftVals: Trying to set deltas gammas for factor " + (*itfac)->getHistName() + "Index out of bounds"); 
			delta = ((prc[j] + prc[j +ll]) * 0.5 - base) * scale_delta;
			gamma = (prc[j] - prc[j + ll]) * scale_gamma;
			del.push_back(delta);
			gam.push_back(gamma);
			++j;
		}
		(*itfac)->setDelta(del);
		(*itfac)->setGamma(gam);
		_pos.push_back(*itfac);
	}
	return "";
}
*/
/*
void readVals_SetSens(string dir_err, string dir_shifts, string dir_pfl_vals, string dir_output, std::vector<boost::shared_ptr<Factor> >& factors, std::vector<boost::shared_ptr<Factor> >& f_Eqvol,
	std::vector<boost::shared_ptr<Factor> >&  f_Ircrv, std::vector<boost::shared_ptr<Factor> >& f_Fxfwd, FactSet& pfl, double shift) {
	
	std::vector<unsigned int> nterms;
	std::vector<std::string> fact_names;
	ErrorLog err_log(dir_err + "err_readVals_SetSens.csv");
	string result, res2;

	// Info files names for Eq_vol, IR-Curve, FX_Forward
	string fn_delta_info("scen_delta_info.csv"), fn_cross_gammas_eqvol_info("scen_cross_gammas_eqvol_info.csv"), fn_cross_gammas_ircrv_info("scen_cross_gammas_ircrv_info.csv"), 
		fn_cross_gammas_fxfwd_info("scen_cross_gammas_fxfwd_info.csv");

	// FactSet valuations for deltas and cross gammas
	std::string fn("pfl_deltas.csv"), fn_gammas_eqvol("pfl_cross_gammas_eqvol.csv"), fn_gammas_ircrv("pfl_cross_gammas_ircrv.csv"), fn_gammas_fxfwd("pfl_cross_gammas_fxfwd.csv");
	std::map<std::string, std::vector<double> > pfl_vals_deltas, pfl_vals_gammas;
		
	// Single term shifts: read valuations, factors info and set sensitivities
	result = readShiftScen_Values(dir_pfl_vals + fn, pfl_vals_deltas);
	err_log.add_line(result);
	result = pfl.setDeltasGammasFromShiftVals(factors, pfl_vals_deltas, shift);
	err_log.add_line(result);

	// Pairs shifts
	// Eq_vol
	pfl_vals_gammas.clear();
	result = readShiftScen_Values(dir_pfl_vals + fn_gammas_eqvol, pfl_vals_gammas);
	err_log.add_line(result);
	result = pfl.setCrossGammasFromShiftVals(f_Eqvol, pfl_vals_gammas, shift);
	err_log.add_line(result);
	
	// Ir_curve
	pfl_vals_gammas.clear();
	result = readShiftScen_Values(dir_pfl_vals + fn_gammas_ircrv, pfl_vals_gammas);
	err_log.add_line(result);
	result = pfl.setCrossGammasFromShiftVals(f_Ircrv, pfl_vals_gammas, shift);
	err_log.add_line(result);

	// Fx_forward
	pfl_vals_gammas.clear();
	result = readShiftScen_Values(dir_pfl_vals + fn_gammas_fxfwd, pfl_vals_gammas);
	err_log.add_line(result);
	result = pfl.setCrossGammasFromShiftVals(f_Fxfwd, pfl_vals_gammas, shift);
	err_log.add_line(result);

	pfl.print_sensitivities(dir_output + "sensitivities.csv");
	err_log.print();
}
*/
/*
std::string FactSet::setCrossGammasFromShiftVals(std::vector<boost::shared_ptr<Factor> >& factors, std::map<std::string, std::vector<double> >& pfl_vals, double shift) {
	std::map<std::string, std::vector<double> >::iterator itmap(pfl_vals.find(_name));

	if(itmap == pfl_vals.end())
		return false;

	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::vector<double> prc(itmap->second);
	int ll((prc.size() - 1) / 4);
	unsigned int iterm, jterm, j(0), nterms;
	double gamma, base(prc[prc.size() - 1]), scale_gamma(1.0 / shift / shift);
	if(ll < 1 || (4 *ll + 1) != prc.size())
		return false;

	for(itfac = factors.begin(); itfac != factors.end(); ++itfac) {
		nterms = (*itfac)->getNterms();
		matrix<double> gam(nterms, nterms);
		for(iterm = 0; iterm < nterms; ++iterm) {
			for(jterm = iterm + 1; jterm < nterms; ++jterm) {
				if(j + 3 * ll >= prc.size()) 
					return std::string("setCrossGammasFromShiftVals: Trying to set cross gammas for factor " + (*itfac)->getHistName() + "Index out of bounds"); 
				gamma = (prc[j] - prc[j + ll] - prc[j + 2 *ll] + prc[j + 3 * ll]) * scale_gamma;
				gam(iterm, jterm) = gamma;
				gam(jterm, iterm) = gamma;
			}
			++j;
		}
		(*itfac)->setCrossGamma(gam);
	}
	return "";
}
*/
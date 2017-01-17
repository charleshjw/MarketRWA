#include "stdafx.h"

//extern std::vector<std::string> factTypeNames;
std::string strNames[] = {"UnIdentified", "FX_Forward", "FX_Spot", "FX_Pair_Vol", "EQ_Vol", "factor_GB", "IRCap_Vol", 
	"IRSwaption_Vol", "CDS", "Turnover", "Refi", "MBS", "IR_Curve", "Eqt_Spot"};
std::vector<std::string> factTypeNames(std::vector<std::string>(strNames, strNames + sizeof(strNames)/sizeof(std::string)));

std::map<std::string, std::string> IrCapsFactorsMap = map_list_of 
	("I07IRVol", ":index_cap1m11") ("I08IRVol", ":index_cap1m21") ("I09IRVol", ":index_cap3m11") ("I10IRVol", ":index_cap3m21") ("I11IRVol", ":index_cap6m11") ("I12IRVol", ":index_cap6m21") ("I02IRVol", ":index_euro1");
std::map<std::string, std::string> FxVolFactorsMap=map_list_of 
	("V03Vol", ":AUD_USD-Simulation") ("V02Vol", ":AUD_JPY-Simulation") ("V22Vol", ":CAD_GBP-Simulation") ("V04Vol", ":CAD_JPY-Simulation") ("V44Vol", ":CAD_USD-Simulation") 
	("V23Vol", ":CHF_GBP-Simulation") ("V08Vol", ":CHF_JPY-Simulation") ("V45Vol", ":CHF_USD-Simulation") ("V25Vol", ":DKK_GBP-Simulation") ("V47Vol", ":DKK_USD-Simulation") 
	("V29Vol", ":GBP_JPY-Simulation") ("V31Vol", ":GBP_NOK-Simulation") ("V33Vol", ":GBP_USD-Simulation") ("V52Vol", ":HKD_USD-Simulation") ("V55Vol", ":JPY_USD-Simulation") 
	("V56Vol", ":MXN_USD-Simulation") ("V58Vol", ":NOK_USD-Simulation") ("V39Vol", ":NZD_USD-Simulation") ("V60Vol", ":RUB_USD-Simulation") ("V62Vol", ":SEK_USD-Simulation") 
	("V63Vol", ":USD_ZAR-Simulation") ("V09Vol", ":CHF_SEK-Simulation") ("V40Vol", ":SEK_DKK-Simulation") ("V73Vol", ":EUR_AUD-Simulation") ("V71Vol", ":EUR_CHF-Simulation") 
	("V67Vol", ":EUR_GBP-Simulation") ("V68Vol", ":EUR_JPY-Simulation") ("V69Vol", ":EUR_NOK-Simulation") ("V70Vol", ":EUR_SEK-Simulation") ("V66Vol", ":EUR_USD-Simulation") 
	("V72Vol", ":EUR_CZK-Simulation") ("V74Vol", ":EUR_TRY-Simulation") ("V75Vol", ":USD_TRY-Simulation");

using boost::math::normal; // typedef provides default type is double.

std::string str_names[] = {"Not Mapped", "Idiosyncratic", "Regressed", "Proxied", "PCA", "NoHistData_NoProxy", "NoHistData_ProxyOfDiffAssetClass","Unknown"};
std::vector<std::string> mapping_names=std::vector<std::string>(str_names, str_names + sizeof(str_names)/sizeof(std::string));

Factor::Factor(std::string s_name) : _nTerms(1), _curr("USD"), _s_name(s_name), _resid_vol(0.0), _R2(0.0), _theTermForMapping(0),_var(0.0),_vol_from_mappedto_facs(0.0) {
	 _h_name = s_name_2_h_name(s_name); 
	 _terms.resize(_nTerms);
}

void Factor::setAttributes(std::map<std::string, std::string>& cur_2_reg_map, std::string hr_data_type, std::string factorGroup, std::string fullName, std::string curr, std::vector<date>& hr_dates, 
	std::map<date, std::vector<double> >& histDataMap, std::map<date, std::vector<unsigned int> >& histTermsMap, int importance) {
		std::string defaultCur("USD"); // Tanya: put it into config file
	unsigned int len(0);
	std::map<date, std::vector<unsigned int> >::iterator itmap;
	std::map<date, std::vector<double> >::iterator itm;
	std::map<std::string, std::string>::iterator it_cur_reg_map;

	FactType iType(UnIdentified);
	_h_name = fullName;
	_hrDates = hr_dates;
	_curr = curr;
	if(_curr.length() > 3)
		_curr = _curr.substr(0, 3);
	_importance = importance;
	_resid_vol = 0.0;
	_R2 = 0.0;
	_vol_from_mappedto_facs = 0.0;
	_var = 0.0;
	_isLogNormal = true;
	_proxy_name = "";
	_histDataAvalType = Unidentified;
	_mapping_type = Unknown;
	_useProxyEOD = false;
	_excludeFromScenarios = false;
	_theTermForMapping = 0;

	// Assign geographical region
	setRegBasedOnCur(cur_2_reg_map);
	
	// Tanya: figure out Factor Type. Hate to do this based on naming conventions
	_type = UnIdentified;
	if(hr_data_type == "hr_spot") {
		if(factorGroup == "fx_data")   {
			if(fullName.find("Spot") != std::string::npos)
				iType = FX_Spot;
			else if(fullName.find("Vol") != std::string::npos)
					if(fullName.find("IRVol") != std::string::npos)
						iType = IRCap_Vol; // or swaption vol?
					else
						iType = FX_Pair_Vol;
		}
		else if(factorGroup ==  "equity_data") {
			if(fullName.find("IRVol") != std::string::npos)
				iType = factor_GB;
		}
		else if(factorGroup ==  "EQUITY_DATA") {
			if(fullName.find("Spot") != std::string::npos)
				iType = Eqt_Spot;
		}
		else if(factorGroup ==  "IR_DATA" && fullName == "USDMBS_CC_OAS") {
			iType = MBS;
			_isLogNormal = false;
		}
		else if(fullName == "REFI_DAILY" || fullName == "TURNOVER_DAILY")
			iType = MBS;
	}

	if(hr_data_type == "hr_curve") {
		if(fullName.substr(0,4) == "CDS_" && fullName.substr(fullName.length()-5, fullName.length()) == "_EQTV")
			iType = CDS;
		else{
			if(factorGroup == "ir_data")   {
				if(fullName.find("Forward") != std::string::npos)
					iType = FX_Forward;
				else {
					iType = IR_Curve;
					_subType = fullName.substr(3, fullName.length());
				}
			}
			if(factorGroup == "EQUITY_DATA")
				if(fullName.find("_EQTV") != std::string::npos)
					iType = EQ_Vol;
		}
	}
	_type = iType;
	_s_name = h_name_2_s_name_by_type(_type);
	_theTermForMapping = 0;
	_hrDataMap = histDataMap;
	_hrTermsMap = histTermsMap;

	len = hr_dates.size();
	if(len > 0) {
		itmap = histTermsMap.find(hr_dates[len-1]);
		_terms = itmap->second;
		_nTerms = _terms.size();
		setEODVal();
	}
	else
		_nTerms = 0;
	
	histDataMap.clear();
	histTermsMap.clear();
	hr_dates.clear();
}
date Factor::setEODVal() {
	bool isFound(false);
	int i;
	unsigned int j;
	date dt(1900,1,1);
	std::vector<unsigned int> terms;
	std::vector<double> eod;
	std::map<date, std::vector<double> >::iterator it;
	double unitHistorical2Eod(1.0); 
	if(_type == FX_Pair_Vol )
		unitHistorical2Eod = 0.01;
	/*
	if(_type == FX_Pair_Vol ) { // Temp Tanya this is how it currently works
		for(i = 0; i < _nTerms; ++i)
			_eod_val.push_back(1.0);
	}
	else if( _hrDates.size() > 0) {
	*/
	if( _hrDates.size() > 0) {
		for(i = _hrDates.size() - 1; i >=0; --i) {
			dt = _hrDates[i];
			it = _hrDataMap.find(dt);
			if(it != _hrDataMap.end()) {
				eod = it->second;
				for(j = 0; j < eod.size(); ++j)
					_eod_val.push_back(eod[j] * unitHistorical2Eod); 
				break;
			}
		}
	}
	else // Temp Tanya: if there is no historical data, set eod to 0 for now. EOD should be read from an eod data file
		for(i = 0; i < _nTerms; ++i)
			_eod_val.push_back(0.0);

	return dt;
}

void Factor::factPrint(std::ofstream& of) {
		of << "Scen name: " << _s_name << "," << "Hist name: " << _h_name << factTypeNames[_type] << "," <<  _curr  << "," << _subType;
		for(unsigned int i = 0; i < _terms.size(); ++i)
			of << "," << _terms[i];

		of << std::endl;
}
void Factor::print_factor_info(std::ofstream& ff) {
	unsigned int iterm, i;
	// Factor name, Mapping type, Factor type, Currency, Region, Number of dates, Number of terms, vol, Resid vol, R2, Regr coef
	ff << _h_name << ", " << mapping_names[_mapping_type] << "," << _type << ", " << _curr << ", " << _reg << ", " << _hrDates.size() << "," << _nTerms << ", " ;
	
	for(iterm = 0; iterm < _nTerms; ++iterm) 
		ff << _vol[iterm]  << ",";
	for(iterm = 0; iterm < 1; ++iterm)
		ff << _R2  << ",";
	for(iterm = 0; iterm < 1; ++iterm)
		ff << _resid_vol   << ",";
	for(i = 0; i < _regr_coef.size1(); ++i)
		ff << _regr_coef(i,0) << ",";
	
	ff << std::endl;
}
void Factor::print_vars(std::stringstream& var_out) {
	unsigned int i;
	double ratio(0.0);
	if(_mapping_type == Idio)
		ratio = 1.0;
	if(_mapping_type == Regressed)
		ratio = pow(_resid_vol/_vol_from_mappedto_facs,2);

	for(i = 0; i < _delta.size(); ++i) {
		var_out << i << "," << _s_name << "," << mapping_names[_mapping_type] << ","  << factTypeNames[_type] << "," << _curr << "," << _reg << ",";
		if(i < _eod_val.size())
			var_out << _eod_val[i];
		var_out << "," << _vol[i] << "," << ratio << "," <<   _delta[i] << ","  << _gamma(i,i) << "," << _var_by_term[i] << "," << _var << std::endl;
	}
}
void Factor::setSpot(const std::vector<double>& prc) {
	unsigned int iterm;
	_spot.resize(prc.size());
	for(iterm = 0; iterm < prc.size(); ++iterm)
		_spot[iterm] = prc[iterm];
}
void Factor::setForecast() {
	// scale historical differences according to forecasted levels and historical volatilities implied from forecast
	// Statistic analysis is needed to forecast FX_Pair_Vol, IRCap_Vol, IRSwaption_Vol, MBS, Refi, Turnover
	unsigned int i, j;
	//std::vector<double> default_scaler(_terms.size(), 1.0);
	//std::vector<double> scaleSpot(_terms.size(), 1.0);
	double scaleSpot(1.0), scaleVol(1.0), tmp;
	switch (_type) {
		case IR_Curve:
			// Need to take care of sub-types that are represented with spreads
		case FX_Forward:
			if(_curr == "USD") {
				// Use the ratio of USDLibor forecast and current USDLibor
			}
			else if(_curr == "GBP") {
				// Use the ratio of GBP Libor forecast and current GBP Libor
			}
			else { // Use the ratio of EUR Libor forecast and current EUR Libor

			}
			break;
		case factor_GB:
		case EQ_Vol:
			// Use the ratio of VIX forecast and current VIX
			break;
		case FX_Spot:
			if(_curr == "USD") {
				// Use the ratio of USDLibor forecast and current USDLibor
			}
			else if(_curr == "GBP") {
				// Use the ratio of GBP Libor forecast and current GBP Libor
			}
			else { // Use the ratio of EUR Libor forecast and current EUR Libor

			}
			break;
		case FX_Pair_Vol:

			break;
		case IRCap_Vol:

			break;
		case IRSwaption_Vol:

			break;
		case MBS:

			break;
		case Refi:

			break;
		case Turnover:

			break;
		//default:
	}
	tmp = scaleSpot * scaleVol;
	for(j = 0; j < _rets.size1(); ++j)
		for(i = 0; i < _terms.size(); ++i)
			_rets(j, i) *= tmp;
	/*
	matrix<double>* ptr = _rets.get();
	for(j = 0; j < (*ptr).size1(); ++j)
		for(i = 0; i < _terms.size(); ++i)
			(*ptr)(j, i) *= tmp;
	*/
}
void Factor::setDelta(std::vector<double>& vals) {
	unsigned int i, n(vals.size());
	_delta.resize(n);

	for(i = 0; i < n; ++i)
		_delta[i] = vals[i];
	//	_delta[i] = 1.0 / _eod_val[i];
}
void  Factor::setDelta(unsigned int iterm, double val) {
	_delta.resize(_nTerms);
	if(iterm < _nTerms)
		_delta[iterm] = val;
}
void Factor::setGamma(std::vector<double>& vals) {
	unsigned int i, n(vals.size());
	zero_matrix<double> zm(n, n);

	_gamma = zm;
	for(i = 0; i < n; ++i)
		_gamma(i, i) = vals[i];
}
void Factor::setCrossGamma(matrix<double>& gam) {
	if(gam.size1() != _gamma.size1() || gam.size2() != _gamma.size2())
		return;
	
	unsigned int n(_gamma.size1()), i, j;

	for(i = 0; i < n; ++i)
		for(j = i + 1; j < n; ++j) {
			_gamma(i, j) = gam(i, j);
			_gamma(j, i) = gam(i, j);
		}
}
/*
boost::shared_ptr<matrix<double> > Factor::getHistData() {
	return _histData;
}

boost::shared_ptr<matrix<double> > Factor::getRets() {
	return _rets;
}
*/
std::string	Factor::h_name_2_s_name() {
	// This is convoluted way to do it. At this point it is hard to do better. Needs to be redone Tanya
	int i;
	std::string s_name(""), h_name(""), curr(""), subType("");
	_s_name = "";
	if(_h_name == "C00891180_TSE_AC___BSpot")
		int tt = 0;

	for(i = 0; i < num_FactType; ++i) {
		s_name = h_name_2_s_name_by_type(i);
		h_name = s_name_2_h_name(s_name);
		if(h_name == _h_name) {
			_s_name = s_name;
			_type = (FactType)i;
			break;
		}
	}
	if(_s_name == "") { // This should happen for either IR_Curve, FX_Forward or FX_Spot
		curr = _h_name.substr(0,3);
		if(_h_name.length() >= 4 && _h_name.substr(_h_name.length() -4, _h_name.length()) == "Spot") {// FXSpot
			_s_name = ":FX" + curr;
			_type = FX_Spot;
		}
		else if(_h_name.length() >= 7 &&_h_name.substr(_h_name.length() - 7, _h_name.length()) == "Forward") {// FX_Forward
			_s_name = ":FX" + curr + "_Forward";
			_type = FX_Forward;
		}
		else if(_h_name.length() >= 3) { // IR_Curve
			subType = _h_name.substr(3, _h_name.length()); 
			if(subType == "MM")
				subType = "Mmarket";
			_s_name = ":IR" + curr + "_" + subType; 
			_type = IR_Curve;
		}
	}
 	return _s_name;
}
std::string  Factor::h_name_2_s_name_by_type(int type) {
	std::map<std::string, std::string>::iterator it;
	std::string s_name(""), subType("");

	switch (type) {
		case IR_Curve: // example: sens ,:IRAUD_Interbank, hist AUDInterbank
			subType = _subType;
			if(subType == "MM")
				subType = "Mmarket";
			s_name = ":IR" + _curr + "_" + subType;
			break;
		case FX_Forward: // Example: sens :FXCAD_Forward, hist CADForward
			s_name = ":FX" + _curr + "_Forward";
			break;
		case FX_Spot: // Example: sens :FXCAD, hist CADSpot
			s_name = ":FX" + _curr;
			break;
		case FX_Pair_Vol: // Example: sens:	:EUR_USD-Simulation, hist V66Vol
			it = FxVolFactorsMap.find(_h_name);
			if(it != FxVolFactorsMap.end())
				s_name = it->second;
			break;
		case EQ_Vol:
			// Example: sens :ACS_C00819010_EQTV, hist ACS_C00819010_EQTV
			s_name = ":" + _h_name;
			break;
		case factor_GB:

			break;
		case Eqt_Spot: // in scenari file "C06738C77_NYSE_DJP", h_name = "C06738C77_NYSE_DJPSpot"
			s_name = ":" + _h_name.substr(0, _h_name.length()-4) + "_Return";
			// This is because of inconsistent naming convention Temp Tanya same logic is in read_prev_day_RF_names
			if(_h_name == "C85299600_OTC_SPTRSpot")
				s_name = ":C85299600_Index_SPTR_Return";
			if(_h_name == "C12496G10_OTC_RUSSELL3000Spot")
				s_name = ":C12496G10_Index_RUSSELL3000_Return";
			if(_h_name == "C00000117_Index_SPALNSpot")
				s_name = ":C00000117_Index_SP500_Return";
			if(_h_name == "C12497K10_OTC_VIXSpot")
				s_name = ":C12497K10_Index_VIX_Return";
			if(_h_name == "C26099405_OTC_INDUASpot")
				s_name = ":C26099405_Index_INDUA_Return";
			break;
		case IRCap_Vol:
		case IRSwaption_Vol:
			it = IrCapsFactorsMap.find(_h_name);
			if(it != IrCapsFactorsMap.end())
				s_name = it->second;
			break;
		
		case CDS: // CDS_RADIOSHA_7C547B_SNRFOR_MR_EQTV
			s_name = ":" + _h_name;
			break;
		case Turnover:
			s_name = ":TURNOVER_DAILY";
			break;
		case Refi:
			s_name = ":Refi_DAILY";
			break;
		case MBS:
			s_name = ":IRUSD_MBS_CC_OAS";
			break;
		//default:
	}

	return s_name;
}

int Factor::comp_Scen_Hist_names(std::string curr, std::string subType, FactType type, std::string h_name) {
	int cmp(1);
	if(_type == type) {  
		switch (_type) {
			case IR_Curve: // example: sens ,:IRAUD_Interbank, hist AUDInterbank
				if(curr == _curr && subType == _subType) {
					_h_name = h_name;
					cmp = 0;
				}
				break;
			case FX_Forward:
			case FX_Spot:
				// Example: sens :FXCAD_Forward, hist CADForward
				if(curr == _curr) {
					_h_name = h_name;
					cmp = 0;
				}
				break;
			
			case FX_Pair_Vol:
				// Example: sens:	, hist V66Vol
				{
				std::map<std::string, std::string>::iterator it = FxVolFactorsMap.find(_s_name);  // todo: change Tanya
				if(it != FxVolFactorsMap.end())
					if(it->second == h_name) {
						_h_name = h_name;
						cmp = 0;
					}
				}
				break;
			case EQ_Vol:
				// Example: sens :ACS_C00819010_EQTV, hist ACS_C00819010_EQTV
				cmp = 0;
				_curr = curr; // Currency cannot be determined from the equity names in scenario file. Currency is specified in historical data file
				break;
			case factor_GB:

				break;
			case IRCap_Vol:

				break;
			case IRSwaption_Vol:

				break;
			case Turnover:

				break;
			case Refi:

				break;
			case MBS:

				break;
			//default:

		}
	}
	
	return cmp;
}

bool Factor::calcPnL_HR(std::vector<double>& pnl) {
	bool isCalced(false);
	if(!_delta.empty()) {
		//if(_spot.empty())
			//file_err_log << "calcPnL_HR: " << "Spot is not available for " << _s_name << std::endl;
		//else if(!_rets)
		//	file_err_log << "calcPnL: " << "Returns are not available for " << _s_name << std::endl;
		if(!_spot.empty()) {
			int i, ndates((_rets).size1());
			for(i = 0; i < ndates; ++i)
				pnl[i] = (_rets)(i, 0) * _spot[0] * _delta[0];
			
			isCalced = true;
		}
	}
	return isCalced;
}
bool Factor::calcPnL_MC(std::vector<double>& vec, unsigned int m0) {
	bool isCalced(false), useGamma(true);
	unsigned int i, j, isim, nsim((_sim_rets).size1()), nterms(_sim_rets.size2());
	boost::numeric::ublas::vector<double> pnl(nsim);
	matrix<double> chg(nsim, nterms);
	double vv(0);
	
	if(_gamma.size1() == nterms && _gamma.size2() == nterms)
		useGamma = true;

	_var_by_term.clear();
	_var_by_term.resize(nterms, 0); 
	vec.clear();
	vec.resize(nsim, 0.0);

	if(_h_name == "USDINAAASP")
		int tt = 0;

	if(_delta.size() == nterms && _eod_val.size() == nterms) {
		for(isim = 0; isim < nsim; ++isim)
			for(i = 0; i < nterms; ++i) {
				// chg(isim, i)  = scaler * _eod_val[i] * (exp(_sim_rets(isim, i) ) - 1.0); 
				// For easier testing temp Tanya
				chg(isim, i)  = _eod_val[i] * _sim_rets(isim, i);
				vec[isim] += _delta[i] * chg(isim, i);
				useGamma = false;
				if(useGamma) {
					for(j = 0; j < nterms; ++j) {
						chg(isim, j)  = _eod_val[j] * _sim_rets(isim, j);
						vec[isim] += 0.5 * _gamma(i,j) * chg(isim,i) * chg(isim, j);
					}
				}
				
			}
		// To calculate VaR by terms for testing mostly
		for(i = 0; i < nterms; ++i) {
			for(isim = 0; isim < nsim; ++isim) {
				chg(isim, i)  = _eod_val[i] * _sim_rets(isim, i);
				pnl[isim] = _delta[i] * chg(isim, i);
				if(useGamma)
					for(j = 0; j < nterms; ++j)
						pnl[isim] += 0.5 * chg(isim, i) *chg(isim, j) * _gamma(i,j);
			}
			sort(pnl.begin(), pnl.end());
			_var_by_term[i] = pnl[m0];
		}

		isCalced = true;
	}
	return isCalced;
}

bool Factor::setRetsVols(ParamControl& params, std::vector<date>& dts) {
	// Finds common dates between historical dates for this factor and input dates, modifies _hrDates, and _hrDataMap, calculates returns and vols
	unsigned int idate, iterm, nrets, minRetsForRegr(params.get_minNumRetsForRegr()), sim_horizon(params.getSimHorizon());
	std::vector<date> dts_cmn(max(dts.size(), _hrDates.size()));
	std::vector<date>::iterator it;
	std::map<date, std::vector<double> >::iterator itmap;
	std::map<date, std::vector<double> > hrDataMapCp(_hrDataMap);
	std::vector<double> prc, ret, wts;
	double vv(1.0), defaultIdioVol(params.get_defaultIdioVol()), defaultIdioEquityVol(params.get_defaultIdioEquityVol()), decay(params.get_decay_factor());

	if(_h_name == "USDMUAAASP")
		int tt = 0; // Temp Tanya
	
	_vol.clear();
	if(_nTerms == 0) {
		if(_type == Eqt_Spot)
			_vol.push_back(defaultIdioEquityVol);
		else
			_vol.push_back(defaultIdioVol);
		_mapping_type = Idio;
		return false;
	}
	// Calc returns and modify dates for this
	setRetsOnDates(dts, sim_horizon, _rets, dts_cmn);
	_hrDates.clear();
	_hrDates = dts_cmn;

	// Set historical data based on modified _hrDates
	_histData.clear();
	_histData.resize(_hrDates.size(), _nTerms);
	_hrDataMap.clear();
	idate = 0;
	for(it = _hrDates.begin(); it != _hrDates.end(); ++it) {
		itmap =hrDataMapCp.find(*it);
		if(itmap != hrDataMapCp.end()) {
			prc = (*itmap).second;
			_hrDataMap.insert(std::pair<date, std::vector<double> >(*it, prc));

			for(iterm = 0; iterm < _nTerms; ++iterm)
				_histData(idate, iterm) = prc[iterm];
	
			++idate;
		}
	}

	// Calc vols
	nrets = _rets.size1();
	if(nrets < minRetsForRegr) {
		for(iterm = 0; iterm < _nTerms; ++iterm) 
			if(_type == Eqt_Spot)
				_vol.push_back(defaultIdioEquityVol);
			else
				_vol.push_back(defaultIdioVol);
		_mapping_type = Idio;
		return false;
	}
	
	ret.resize(nrets);
	calculate_weights(decay, _hrDates, wts, sim_horizon);
		
	for(iterm = 0; iterm < _nTerms; ++iterm) {
		for(idate = 0; idate < nrets; ++idate) 
			ret[idate] = _rets(idate, iterm);
		
		_vol.push_back(stdev(ret, wts));
	}
	return true;
}
bool Factor::setRetsOnDates(std::vector<date>& dts, unsigned int sim_horizon, matrix<double>& rets, std::vector<date>& dts_cmn) {
	// Calculates returns on dates and returns the result without modifying neither _rets, nor _hrDates and _hrDataMap
	std::vector<date>::iterator it;
	std::map<date, std::vector<double> >::iterator itmap;
	std::vector<double> vec;
	unsigned int iterm, idate(0);
	date date1, date2;
	date_duration dt(0);
	matrix<double> prc;
	unsigned int nrets;
	double t1, v1, eps(1E-7);

	rets.clear();
	dts_cmn.resize(max(dts.size(), _hrDates.size()));
	it = set_intersection(dts.begin(), dts.end(), _hrDates.begin(), _hrDates.end(), dts_cmn.begin());
	dts_cmn.resize(it - dts_cmn.begin());

	if(dts_cmn.size() <= sim_horizon + 1 || _nTerms == 0)
		return false;

	nrets = dts_cmn.size() - sim_horizon;

	prc.resize(dts_cmn.size(), _nTerms);
	rets.resize(nrets, _nTerms);
	
	for(idate = 0; idate < dts_cmn.size(); ++idate) {
		itmap =_hrDataMap.find(dts_cmn[idate]);
		vec = (*itmap).second;
		for(iterm = 0; iterm < _nTerms; ++iterm)
			prc(idate, iterm) = vec[iterm];
	}

	for(idate = sim_horizon; idate < dts_cmn.size(); ++idate) {
		dt = dts_cmn[idate] - dts_cmn[idate - sim_horizon];
		for(iterm = 0; iterm < _nTerms; ++iterm)
			if(_isLogNormal) {
				v1 = 0.0;
				if(abs(prc(idate - sim_horizon, iterm)) > eps) {
					t1 =  prc(idate, iterm) / prc(idate - sim_horizon, iterm);
					if(t1 > 0)
						v1 = log(t1) / sqrt(dt.days() * 1.0);
				}
				rets(idate - sim_horizon, iterm) = v1;
			}
			else
				rets(idate - sim_horizon, iterm) = (prc(idate, iterm) - prc(idate - sim_horizon, iterm)) / sqrt(dt.days() * 1.0);
	}

	return true;
}

bool Factor::mapIt(ParamControl& params, std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Spot, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Spot, 
	std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_IR_Curve,
	std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Forward, 	std::vector<date>& dts, std::map<std::string, std::vector<std::string> >& reg_2_cur_map) {
	
	double decay(params.get_decay_factor()), defaultIdioVol(params.get_defaultIdioVol()), defaultIdioEquityVol(params.get_defaultIdioEquityVol());
	bool success(false), enoghData2Regr(false);
	std::vector<boost::shared_ptr<Factor> > factors_4_regression;
	std::vector<std::string> currencies_4_mapping;
	std::vector<date>::iterator itdt;
	unsigned int sim_horizon(params.getSimHorizon());

	_mapping_type = Idio;
	_vol.clear();
	if(_type == Eqt_Spot)
		_vol.push_back(defaultIdioEquityVol);
	else
		_vol.push_back(defaultIdioVol);

	if(_h_name == "C00287Y10_NYSE_ABBVSpot" )
		int tt = 0;

	enoghData2Regr = setRetsVols(params, dts);

	if(!enoghData2Regr) 
		return enoghData2Regr;

	// Find set of factors to regress against
	switch(_type) {
		case IR_Curve:	// Regress against data for the same region
			/*
			for(itcur = currencies_4_mapping.begin(); itcur != currencies_4_mapping.end(); ++itcur) {
				jt = findFactor(_type, _subType, *itcur, f_IR_Curve);
				if(jt != f_IR_Curve.end())
					factors_4_regression.push_back((*jt));
			}
			*/
			factors_4_regression = f_PCA_IR_Curve;
			break;
		case FX_Forward:
			/*
			for(itcur = currencies_4_mapping.begin(); itcur != currencies_4_mapping.end(); ++itcur) {
				jt = findFactor(_type, _subType, *itcur, f_FX_Forward);
				if(jt != f_FX_Forward.end())
					factors_4_regression.push_back((*jt));
			}
			*/
			factors_4_regression = f_PCA_FX_Forward;
			break;

		case FX_Spot:
			factors_4_regression = f_PCA_FX_Spot;
			/*
			for(itcur = currencies_4_mapping.begin(); itcur != currencies_4_mapping.end(); ++itcur) {
				jt = findFactor(_type, _subType, *itcur, f_FX_Spot);
				if(jt != f_FX_Spot.end())
					factors_4_regression.push_back((*jt));
			}
			*/
			break;
		case Eqt_Spot: // Regress agains indices
		case EQ_Vol:
			factors_4_regression = f_PCA_Eqt_Spot;
			break;

		case CDS:
			break;
		case FX_Pair_Vol:
			factors_4_regression = f_PCA_FX_Vol;
			break;
		
		case IRCap_Vol:
		case IRSwaption_Vol:
		default:
			break;
		// Turnover, Refi, MBS, factor_GB, 
	}
	if(factors_4_regression.size() == 0) 
		return false;
	
	// Regression
	success = regress(factors_4_regression, sim_horizon, decay);

	return success;
}
bool Factor::regress(std::vector<boost::shared_ptr<Factor> >& fvec, unsigned int sim_horizon, double decay) {
	bool success(false);
	unsigned int j, i, nrets(_rets.size1()), nfacs(fvec.size()), ifac, theTermIdx, theTerm;
	std::vector<boost::shared_ptr<Factor> >:: iterator jt;
	std::vector<date>::iterator itdt;
	std::vector<unsigned int> terms;
	std::vector<date> new_dts;
	std::map<date, std::vector<double> > hrDataMap;
	matrix<double> X(nrets, nfacs), Y(nrets,1), Y1(nrets,1), resid(nrets,1), rets;
	double std(0.0), variance_from_mapped_facs(0.0);
	date_duration dt(0);
	std::vector<double> wts, vol_on_all_dates(nfacs), vol_on_these_dates(nfacs);

	_mapping_type = Idio;
	if(nrets <= sim_horizon)
		return false;

	theTermIdx = _nTerms * 0.5;
	theTerm = _terms[theTermIdx];
	_theTermForMapping = theTermIdx;

	// Collect historical data for regression
	//_regr_coef.clear();
	//_regr_coef.resize(1, nfacs);
	
	X.resize(nrets, nfacs);
	ifac = 0;
	for(jt = fvec.begin(); jt != fvec.end(); ++jt) {
		// Which term matches theTerm?
		int dd(1000000); // large number
		unsigned int jj(0), nt((*jt)->getNterms());
		terms = (*jt)->getTerms();

		for(i = 0; i < nt; ++i) {
			int vv(terms[i] - theTerm);
			if(abs(vv) < dd) {
				jj = i;
				dd = abs(vv);
			}
		}

		(*jt)->setRetsOnDates(_hrDates, sim_horizon, rets, new_dts);
		if(rets.size1() != Y.size1() || jj >= rets.size2()) {
			return false;
		}
		if(jt == fvec.begin())
			calculate_weights(decay, new_dts, wts, sim_horizon);

		for(j = 0; j < nrets; ++j)
			X(j, ifac) = rets(j, jj) * wts[j];

		_mapped_fact.push_back(*jt);
		_mapped_terms.push_back(jj);
		++ifac;
	}

	
	for(j = 0; j < nrets; ++j)
		Y(j, 0) = _rets(j, theTermIdx) * wts[j];

	// Each column of Y is regressed against X
	unsigned int n(X.size1()), m(X.size2()), k(1);

	if(n == 0 || m == 0 || k == 0)
		return false;

	matrix<double> aa(m, m), bb(m,k), invm(m,m);

	_regr_coef.clear();
	_regr_coef.resize(m, k);
	i = 0;
	aa = prod(trans(X), X);
	for(j = 0; j < nfacs; ++j) {
		vol_on_all_dates[j] = fvec[j]->getVol()[0];
		vol_on_these_dates[j] = sqrt(aa(j, j));
	}

	if(InvertMatrix(aa, invm)) {
		bb = prod(trans(X), Y);
		_regr_coef = prod(invm, bb);
		Y1 = prod(X, _regr_coef);
		resid = Y - Y1;
		
		//_resid_vol = sqrt(inner_product(column(resid,i).begin(), column(resid,i).end(), column(resid, i).begin(), 0.0) / n);
		_resid_vol = sqrt(inner_product(column(resid,i).begin(), column(resid,i).end(), column(resid, i).begin(), 0.0));
		//std = sqrt(inner_product(column(Y,i).begin(), column(Y,i).end(), column(Y,i).begin(), 0.0) / n);
		std = sqrt(inner_product(column(Y,i).begin(), column(Y,i).end(), column(Y,i).begin(), 0.0));
		_R2 = 1.0 - pow(_resid_vol / std, 2);
	}
	variance_from_mapped_facs = inner_product(column(Y1,i).begin(), column(Y1,i).end(), column(Y1,i).begin(), 0.0);
	_vol_from_mappedto_facs = sqrt(variance_from_mapped_facs);
	for(j = 0; j < nfacs; ++j)
		_regr_coef(j, 0) *= (vol_on_these_dates[j]/vol_on_all_dates[j]);
	_mapping_type = Regressed;
	return true;
}
void Factor::setAttributes(FactType iType, std::string subType, std::string cur, std::string reg, std::vector<double>& sens, std::vector<unsigned int> terms, int importance) {
	_type = iType;
	_subType = subType;
	_curr = cur;
	_reg = reg;
	_importance = importance;
	setTerms(terms);
	_delta = sens;
}
void Factor::getTermsAtDate(date dd, std::vector<unsigned int>& terms) {
	
	terms.clear();
	std::map<date, std::vector<unsigned int> >::iterator itmap;
	itmap = _hrTermsMap.find(dd);
	if(itmap != _hrTermsMap.end())
		terms = itmap->second;

	return;
}
void Factor::getDataAtDate(date dd, std::vector<double>& vals) {
	std::map<date, std::vector<double> >::iterator itmap;
	vals.clear();
	itmap = _hrDataMap.find(dd);
	if(itmap != _hrDataMap.end())
		vals = itmap->second;
	
	return;
}
bool Factor::dataCleanUp(ParamControl& params) {
	// Historical data are cleaned up so that for each date we have vector of prices that correspond to the same tenors as tenors available for the last historical date.
	// This is not the robast approach - it is possible that the tenors for the last available date are messed up... Tanya
	// Get read of zero, stale and sparse data
	double eps(1E-6), unreasonableValue(params.get_unreasonableValue()), staleDataDailyVolMin(params.get_staleDataDailyVolMin()), smallValue(0.0001); // smallValue is assigned insead of zero values in historical data files
	unsigned int i, idate(0), ndates(_hrDates.size()), theTerm(0);
	unsigned int iterm, maxTerms(100), minPositiveValue(params.get_minPositiveValue()), maxDatesUnch(params.get_maxDatesUnch()), minDatesCount(params.get_minDatesCount()), maxAllowedGap(params.get_maxAllowedGap());
	std::vector<double> prc;
	std::vector<unsigned int> tr, tr_res;
	std::vector<unsigned int>::iterator iti;
	std::vector<date>::iterator itdt;
	std::vector<date> hrDatesCp(_hrDates);
	std::vector<double> eod;
	double vol;
	bool badData(false);
	matrix<double> data;
	matrix<unsigned int> trms;
	std::map<date, std::vector<double> >::iterator itmap;
	std::map<date, std::vector<unsigned int> >::iterator ittrms;
	double t1, t2;

	_regr_coef.clear();
	_regr_coef.resize(1, 1);

	if(_h_name == "BRKB_C08467070_EQTV" )
		int tt = 0;

	if (ndates == 0)
		return false;
	
	// Read terms data
	idate = 0;
	trms.resize(ndates, maxTerms);
	trms.clear();
	for(ittrms = _hrTermsMap.begin(); ittrms != _hrTermsMap.end(); ++ittrms) {
		tr = ittrms->second;
		if(idate == ndates - 1) {
			_nTerms = tr.size();
			_terms = tr;
		}
		
		for(iterm = 0; iterm < tr.size(); ++iterm) 
			trms(idate, iterm) = tr[iterm];
		
		++idate;
	}

	// Get data into matrix form
	// If the set of terms for some date is different from the set of terms for the last available date then this date is eliminated
	_hrDates.clear();
	data.resize(ndates, maxTerms);
	tr_res.resize(maxTerms);
	idate = 0;
	for(itmap = _hrDataMap.begin(); itmap != _hrDataMap.end(); ++itmap) {
		tr.clear();
		tr.resize(maxTerms);
		i = 0;
		for(iterm = 0; iterm < maxTerms; ++iterm) {
			if(iterm == 0 || (iterm > 0 && trms(idate, iterm) > 0)) {
				tr[i] = (trms(idate, iterm));
				++i;
			}
		}
		tr.resize(i);
		iti = set_intersection(tr.begin(), tr.end(), _terms.begin(), _terms.end(), tr_res.begin());
		tr_res.resize(iti - tr_res.begin());
		prc = itmap->second;
		if(tr_res.size() == _nTerms && prc.size() == _nTerms) {
			for(iterm = 0; iterm < _nTerms; ++iterm)
				data(idate, iterm) = prc[iterm];
			_hrDates.push_back(hrDatesCp[idate]);
			++idate;
		}
	}

	// By now the terms are established
	ndates = idate;
	data.resize(ndates, _nTerms);

	// Identify non-positive of stale data
	bool unch(true);
	for(iterm = 0; iterm < _nTerms; ++iterm) {
		prc.clear();
		for(idate = 0; idate < ndates; ++idate) {
			//if(_isLogNormal && data(idate, iterm) >= 0 && data(idate, iterm) < minPositiveValue)
			//	data(idate, iterm) = minPositiveValue; 
			if(idate >= maxDatesUnch) {
				unch = true;
				for(i = 1; i < maxDatesUnch; ++i) {
					t1 = data(idate - i, iterm);
					t2 = data(idate - i + 1, iterm);
					/*
					if(_isLogNormal && (data(idate - i, iterm) <= minPositiveValue || data(idate - i + 1, iterm) <= minPositiveValue ))
						break;
					
					if(_isLogNormal) 
						vv = log(data(idate - i, iterm)/data(idate - i + 1, iterm));
					else
						vv = data(idate - i, iterm) - data(idate - i + 1, iterm);
					if(abs(vv) > eps) {
						unch = false;
						break;
					}
					*/
					if(abs(t1 - t2) > eps) {
						unch = false;
						break;
					}
				}
				if(unch)
					for(i = 1; i <= maxDatesUnch; ++i)
						data(idate - i, iterm) = unreasonableValue;
			}
			//if(_isLogNormal && data(idate, iterm) <= minPositiveValue)
			//	data(idate, iterm) = unreasonableValue;
			// for test only
			prc.clear();
			for(i=0; i <= idate; ++i)
				prc.push_back(data(i, iterm));
		}
	}
	
	// Select theTerm with the most of good data 
	unsigned int maxCount(0), minCount(100000000), cnt;
	theTerm = 0;
	double ave, ret;
	std::vector<unsigned int> countGoodValues(_nTerms);
	for(iterm = 0; iterm < _nTerms; ++iterm) {
		ave = 0.0;
		vol = 0.0;
		cnt = 0;
		for(idate = 0; idate < ndates; ++idate) {
			if(data(idate, iterm) != unreasonableValue) 
				countGoodValues[iterm] +=1;
			if(idate > 0 ) {
				ret = 0.0;
				if(_isLogNormal) {
					if(abs(data(idate-1, iterm)) > eps) {
						t1 = data(idate, iterm) / data(idate-1, iterm);
						if(t1 > 0.0)
							ret = log(t1);
					}
				}
				else
					ret = data(idate, iterm) - data(idate-1, iterm);
				ave += ret;
				vol += ret * ret;
				++cnt;
			}
		}
		if(cnt > 0) {
			ave /= cnt;
			vol /= cnt;
			vol = sqrt(vol - ave * ave);
		}
		if(cnt < maxDatesUnch || vol < staleDataDailyVolMin) {
			countGoodValues[iterm] = 0;
			for(idate = 0; idate < ndates; ++idate)
				data(idate, iterm) = unreasonableValue;
		}
		
		if(countGoodValues[iterm] > maxCount) {
			maxCount = countGoodValues[iterm];
			theTerm = iterm;
		}
		if(countGoodValues[iterm] < minCount)
			minCount = countGoodValues[iterm];
	}

	// Now  "theTerm" has maximum "good dates" across terms which number is maxCount, minimum "good dates" across tenors is minCount
	matrix<double> dataCp(ndates, _nTerms);
	tr = _terms;
	eod = _eod_val;
	i = 0;
	for(iterm =0; iterm < _nTerms; ++iterm)
		if(countGoodValues[iterm] == maxCount) {
			column(dataCp, i) = column(data, iterm);
			_terms[i] = tr[iterm];
			_eod_val[i] = eod[iterm];
			++i;
		}
	dataCp.resize(ndates, i);
	_terms.resize(i);
	_eod_val.resize(i);
	data = dataCp;
	_nTerms = i;
	/*
	if(minCount < minDatesCount) { // Only theTerm stays
		tr = _terms;
		eod = _eod_val;
		_terms.resize(1);
		_terms[0] = tr[theTerm];
		_nTerms = 1;
		_eod_val.resize(1);
		_eod_val[0] = eod[theTerm];

		for(idate = 0; idate < ndates; ++idate)
			dataCp(idate, 0 ) = data(idate, theTerm);
		data = dataCp;
	}
	*/
	_hrDataMap.clear();
	hrDatesCp = _hrDates;
	_hrDates.clear();
	bool goodDate(true);
	i = 0;
	prc.clear();
	prc.resize(_nTerms);
	dataCp.clear();
	dataCp.resize(data.size1(), data.size2());

	for(idate = 0; idate < ndates; ++idate) {
		goodDate = true;
		for(iterm = 0; iterm < _nTerms; ++iterm) 
			if(data(idate, iterm) == unreasonableValue) {
				goodDate = false;
				break;
			}
		if(goodDate) {
			for(iterm = 0; iterm < _nTerms; ++iterm) 
				prc[iterm] = data(idate, iterm);
		}
		if(goodDate) {
			unsigned int gap(0);
			if(_hrDates.size() > 0) {
				date_duration ll(hrDatesCp[idate] - _hrDates[_hrDates.size() - 1]);
				gap = ll.days();
			}
			if(gap > maxAllowedGap) {
				_hrDates.clear();
				_hrDataMap.clear();
				i = 0;
			}
			_hrDataMap.insert(std::pair<date, std::vector<double> > (hrDatesCp[idate], prc));
			_hrDates.push_back(hrDatesCp[idate]);
			for(iterm = 0; iterm < _nTerms; ++iterm) 
				dataCp(i, iterm) = prc[iterm];
			++i;
		}
	}
	ndates = i;
	dataCp.resize(ndates, _nTerms);

	// Identify large erroneous returns 
	if(ndates > 3) {
		double curChange, nextChange, nextChange2, nextChange3, tol(params.get_maxAllowedRet());
		for(iterm = 0; iterm < _nTerms; ++iterm) {
			matrix_column<matrix<double> > value(dataCp, iterm);
			for(i = 1; i < ndates - 3; ++i) {
				if(_isLogNormal) {
					curChange = value[i] / value[i-1] - 1.0;
					nextChange = value[i + 1] / value[i] - 1.0;
					nextChange2 = value[i + 2] / value[i + 1] - 1.0;
					nextChange3 = value[i + 3] / value[i + 2] - 1.0;
				} else { // For now the only case of not log-normal is USDMBS_CC_OAS and its historical data are the daily differences already
					curChange = value[i] - value[i-1];
					nextChange = value[i + 1] - value[i];
					nextChange2 = value[i + 2] - value[i + 1];
					nextChange3 = value[i + 3] - value[i + 2];
				}
				if ( _isLogNormal && fabs(curChange) > tol ) {		// scub condition is today and tomorrow
					if ( nextChange * curChange < 0.0 				// changes in different direction, while
						&& fabs(nextChange) > tol ) {				// both changes are greater than tolerance
							value[i] = (value[i-1] + value[i+1]) * 0.5;
					} else if ( nextChange2 * curChange < 0.0 && fabs(nextChange2) > tol ) {
						value[i] = (value[i - 1] + value[i+2]) * 0.5;
						value[i+1] = value[i];
					} else if ( nextChange3 * curChange < 0.0 && fabs(nextChange3) > tol ) {
						value[i] = (value[i-1] + value[i+3]) * 0.5;
						value[i+1] = value[i];
						value[i+2] = value[i];
					}
				}
			}
		}
	}
	// Update hist data map
	idate = 0;
	for(itmap = _hrDataMap.begin(); itmap != _hrDataMap.end(); ++itmap) {
		prc.clear();
		for(i = 0; i < _nTerms; ++i)
			prc.push_back(dataCp(idate, i));
		itmap->second = prc;
		++idate;
	}
	return true;
}

void Factor::print_sim_rets(std::ofstream& ff) {
	unsigned int iterm, isim, nsim(_sim_rets.size1()), nterms(_sim_rets.size2());
	std::string factor_type(mapping_names[_mapping_type]);
	
	if(_h_name == "C78462F10_NYSE_SPY_ETFSpot")
		int tt = 0;

	if(nterms == 0) {
		ff << factor_type << "," << _h_name << "," << _s_name  << "," << _subType << "," << factTypeNames[_type] << "," << nterms; 

		if(nterms != _rets.size2())
			ff << "," << _h_name << "," << _rets.size2() << ",Hist rets";

		ff << std::endl;
	}
	for(iterm = 0; iterm < nterms; ++iterm) {
		ff << factor_type << "," << _h_name << "," << _s_name  << "," << _subType << "," << factTypeNames[_type] << "," << nterms << "," << _terms[iterm] << "," ;
		for(isim = 0; isim < nsim; ++isim)
			ff << _sim_rets(isim, iterm) << ",";
		
		ff << std::endl;
	}
}
void Factor::print_sim_diff(std::ofstream& ff) {
	unsigned int iterm, isim, nsim(_sim_rets.size1()), nterms(_sim_rets.size2());
	std::string factor_type(mapping_names[_mapping_type]);
	std::vector<double> vec(nsim);
	double mean(0), vol(1.0);

	if(_h_name == "C78462F10_NYSE_SPY_ETFSpot")
		int tt = 0;

	if(nterms == 0) {
		ff <<  _h_name << "," << _s_name  << "," << _nTerms << "," << _curr << "," << _reg << "," << factor_type << "," << factTypeNames[_type] << "," << _isLogNormal << "," << nterms;

		if(nterms != _rets.size2())
			ff << "," << _h_name << "," << _rets.size2() << ",Hist rets";

		ff << std::endl;
	}
	for(iterm = 0; iterm < nterms; ++iterm) {
		ff << _h_name << "," << _s_name  << "," << _nTerms << "," << _curr << "," << _reg << "," << factor_type << "," << factTypeNames[_type] << "," << _isLogNormal << "," << iterm;
		if(iterm < _terms.size())
			ff << "," << _terms[iterm] << ",";
		else
			ff << ",0,";
		
		ff << _rets.size1() << ",";

		if(iterm < _dist_flags.size() && _dist_flags[iterm])
			ff << "Real Dist,";
		else
			ff << "Normal Dist,";
		
		if(iterm < _vol.size())
			ff << _vol[iterm];
		ff << ",";

		if(iterm < _eod_val.size()) {
			ff << _eod_val[iterm] << ",";
			if(_type == FX_Pair_Vol )
				for(isim = 0; isim < nsim; ++isim)
					vec[isim] = exp(_sim_rets(isim, iterm)) - 1.0;
			else {
				if(_isLogNormal)
					for(isim = 0; isim < nsim; ++isim)
						vec[isim] = _eod_val[iterm] * (exp(_sim_rets(isim, iterm)) - 1.0);
				else
					for(isim = 0; isim < nsim; ++isim)
						vec[isim] = _sim_rets(isim, iterm);
			}
			normalize_vec(vec, mean, vol);
			ff << vol << "," ; 
		}
		else
			ff << ",,";

		for(isim = 0; isim < nsim; ++isim)
			ff << _sim_rets(isim, iterm) << ",";
		
		ff << std::endl;
	}
}
void Factor::print_rets(std::ofstream& ff) {
	unsigned int iterm, isim, nhist(_rets.size1()), nterms(_sim_rets.size2());
	std::string factor_type(mapping_names[_mapping_type]); 

	if(_h_name == "INDUA_C26099405_EQTV")
		int tt = 0;

	if(nterms == 0) {
		ff << "," << factor_type << "," << _h_name << "," << _s_name  << "," << _subType << "," << factTypeNames[_type] << "," << nterms;

		if(nterms != _rets.size2())
			ff << "," << _h_name << "," << _rets.size2() << ",Hist rets";

		ff << std::endl;
	}
	for(iterm = 0; iterm < nterms; ++iterm) {
		if(iterm < _vol.size())
			ff << _vol[iterm];
			
		ff	<< "," << factor_type << "," << _h_name << "," << _s_name  << "," << _subType << ","  << factTypeNames[_type];
		if(iterm < _terms.size())
			ff << "," << _terms[iterm] << ",";
		else
			ff << ",0,";

		ff << nhist << ",";

		if(iterm < _rets.size2())
			for(isim = 0; isim < nhist; ++isim)
				ff << _rets(isim, iterm) << ",";
				
		ff << std::endl;
	}
}

void Factor::print_hist(std::ofstream& ff) {
	unsigned int iterm, ndates(_hrDates.size());
	std::map<date, std::vector<double> >::iterator itmap;
	std::vector<double> prc;
	date dt;
	ff << "Name, date, val" << std::endl;
	ff << _h_name <<"," << ndates << ",,,,";
	for(iterm = 0; iterm < _nTerms; ++iterm)
		ff << _terms[iterm] << ",";
	ff << std::endl;
	
	for(itmap = _hrDataMap.begin(); itmap != _hrDataMap.end(); ++itmap) {
		dt = (*itmap).first;
		ff << _h_name << "," << dt.month() << "-" << dt.day() << "-" << dt.year() << "," << dt.month() << "," <<  dt.day() << "," << dt.year() << ",";
		prc = (*itmap).second;
		for(iterm = 0; iterm < _nTerms; ++iterm)
			ff << prc[iterm] << ",";
		ff << std::endl;
	}
}
void Factor::print_hist_pca(std::ofstream& ff) {
	unsigned int ndates(_hrDates.size());
	std::map<date, std::vector<double> >::iterator itmap;
	std::vector<double> prc;
	date dt;

	ff << _h_name <<"," << ndates << ",";
	for(itmap = _hrDataMap.begin(); itmap != _hrDataMap.end(); ++itmap) {
		dt = (*itmap).first;
		//ff << _h_name << "," << dt.month() << "-" << dt.day() << "-" << dt.year() << ",";
		prc = (*itmap).second;
		ff << prc[0] << ",";
	}
	ff << std::endl;

}
void Factor::calcAdhocShiftTilt(boost::numeric::ublas::vector<double>& shift, boost::numeric::ublas::vector<double>& tilt) {
	// Adhoc formula implies more or less realistic correlations across terms
	unsigned int i;
	const unsigned int n1(356);
	const double c1(0.01);

	shift.resize(_nTerms);
	tilt.resize(_nTerms);
	for(i = 0; i < _nTerms; ++i) {
		shift[i] = exp(-c1 * pow(i, 2.0));
		tilt[i] = sqrt(1.0 - shift[i] * shift[i]);
		if(_terms[i] > n1)
			tilt[i] = -tilt[i];
	}
}
void Factor::genSims(matrix_range<matrix<double> >& rand_norm, std::vector<double>& wts) {
	unsigned int nRets(_rets.size1()), nterms(_rets.size2()), i, j, k, nsim(rand_norm.size1());
	
	matrix<double> weighted(nRets, nterms), sim_rets;
	std::vector<double> empir_vec(nRets), norm_vec(nsim), res(nsim), equal_wts(nsim), sim_cum(nsim), sim_transformed(nsim), empir_x(nRets), empir_y(nRets);
	double percentile(0.01), quantile(2.32635), scaler(1.0), sim_vol(1.0), ww(sqrt(1.0/nsim));
	unsigned int n1(floor(ceil(nRets * percentile))), n2(floor(nRets * (1 - percentile))), m1(ceil(nsim * percentile)), m2(floor(nsim * (1 - percentile)));
	std::vector<unsigned int> norm_vec_idx(nsim), empir_vec_idx(nRets);

	for(i = 0; i < nRets; ++i)
		row(weighted, i) = wts[i] * row(_rets, i);
	
	for(i = 0; i < nsim; ++i)
		equal_wts[i] = ww;

	sim_rets = prod(rand_norm, weighted);

	if(_h_name == "AEDSpot" || _h_name == "EURInterbank")
		int tt = 0;

	_sim_rets = sim_rets;

	// Transform to real dist if needed
	for(i = 0; i < nterms; ++i) {
		for(j = 0; j < nsim; ++j) {
			norm_vec[j] =_sim_rets(j, i);
			sim_transformed[j] = norm_vec[j];
		}

		if(_dist_flags.size() == _nTerms && _dist_flags[i] ) {
			for(j = 0; j < nRets; ++j)
				empir_x[j] = _rets(j, i);
			
			sort(empir_x.begin(), empir_x.end());
			for(k = 0; k < nRets; ++k) 
				empir_y[k] = (k+1.0) / (nRets + 1);
				
			sort(norm_vec.begin(), norm_vec.end());
			for(j = 0; j < nsim; ++j) {
				k = 0;
				while(_sim_rets(j, i) > norm_vec[k])
					++k;
				sim_cum[j] = (k+1.0) / (nsim + 1);
			}
			if(0) { // linear interp
				std::map<double, double> empir_cum;
				for(j = 0; j < nRets; ++j)
					empir_cum.insert(std::make_pair<double, double>(empir_y[j], empir_x[j]));
				interpolate_linear(empir_cum, sim_cum, sim_transformed);
			}
			if(1) { // quadratic
				interpolate_quadratic(empir_y, empir_x, sim_cum, sim_transformed);
			}
			for(j = 0; j < nsim; ++j)
				_sim_rets(j, i) = sim_transformed[j];
			/*
			sort(empir_vec.begin(), empir_vec.end());
			scaler = max(fabs(empir_vec[n1] / norm_vec[m1]), fabs(empir_vec[n2] / norm_vec[m2]));
			column(_sim_rets, i) *= scaler;
			*/
		}
		else {
			sim_vol = stdev(norm_vec, equal_wts);
			scaler = _vol[i]/sim_vol;
			column(_sim_rets, i) *= scaler;
		}
		
	}
}

void Factor::complete(std::vector<date>& dts) {
	unsigned int i, j, iterm, j1, j2;
	double vv(0.0);
	std::vector<date>::iterator it;
	std::vector<double> prc1(_nTerms), prc2(_nTerms), prc(_nTerms);
	std::map<date, std::vector<double> >::iterator itmap;

	if(_h_name == "C06738C77_NYSE_DJPSpot")
		int tt = 0;

	for(i = 0; i < dts.size(); ++i) {
		it = find(_hrDates.begin(), _hrDates.end(), dts[i]);
		if(it == _hrDates.end()) {
			if(dts[i] < _hrDates[0]) {
				j1 = 0;
				j2 = j1;
			}
			else if(dts[i] > _hrDates[_hrDates.size() - 1]) {
				j1 = _hrDates.size() - 1;
				j2 = j1;
			}
			else {
				j = 0;
				while(dts[i] > _hrDates[j])
					++j;
				j1 = j - 1;
				j2 = j;

			}
			itmap = _hrDataMap.find(_hrDates[j1]);
			prc1 = itmap->second;
			if(j1 == j2) 
				prc = prc1;
			else {
				itmap = _hrDataMap.find(_hrDates[j2]);
				prc2 = itmap->second;
			
				date_duration dt1(dts[i] - _hrDates[j1]), dt(_hrDates[j2] - _hrDates[j1]);
				vv = (dt1.days() * 1.0) / dt.days();
				for(iterm = 0; iterm < _nTerms; ++iterm)
					prc[iterm] = prc1[iterm] + vv * (prc2[iterm] - prc1[iterm]);
			}
			_hrDataMap.insert(std::make_pair(dts[i], prc));
			if(dts[i] > _hrDates[_hrDates.size() - 1])
				_hrDates.push_back(dts[i]);
			else
				_hrDates.insert(_hrDates.begin() + j2, 1, dts[i]);
		}
	}
}
void Factor::print_sensitivities(std::stringstream& ff) {
	unsigned int i, j, n1(_delta.size()), n2(_gamma.size1());

	if(n1 != n2)
		return;

	ff << _s_name << ", Delta";
	for(i = 0; i < n1; ++i)
		ff << "," << _delta[i];
	ff << std::endl;

	for(i = 0; i < n2; ++i) {
		ff << _s_name << ", Gamma";
		for(j = 0; j < n2; ++j)
			ff << "," << _gamma(i, j);
		ff << std::endl;
	}
}

double Factor::interpolate(double x, std::string attribute) {
	// Linear interpolation
	double y(0);
	unsigned int i(0), nterms(0);

	std::map<double, double> mp;
	std::map<double, double>::iterator it, jt;
	if(attribute == "EOD") {
		if(_eod_val.size() != _nTerms)
			if(_eod_val.size() > 0)
				return _eod_val[0];
			else
				return 0.0; // Need to handle this Tanya

		for(i = 0; i <_nTerms; ++i)
			mp.insert(std::make_pair(_terms[i], _eod_val[i]));

		it = mp.upper_bound(x);
		if(it == mp.end())
			y = (--it)->second;
		else if(it == mp.begin())
			y = it->second;
		else {
			jt = it;
			--jt;
			y = jt->second + (it->second - jt->second) / (it->first - jt->first) * (x - jt->first);
		}
	}
	return y;
}
void Factor::handle_noData_noProxy(double defaultIdioVol, double defaultIdioEquityVol) {
	_eod_val.resize(1);
	_vol.resize(1);
	_eod_val[0] = 0.0;
	_vol[0] =defaultIdioVol;
	if(_type == Eqt_Spot)
		_vol[0] =defaultIdioEquityVol;
	_nTerms = 1;
	_terms.resize(1);
	_terms[0] = 0;
	_proxy_name ="";
	_mapping_type = Idio;
}
void Factor::scale_eod(std::map<double,double>& map) {

	unsigned int iterm;
	std::vector<double> terms_d(_nTerms), coef(_nTerms);

	if(_eod_saved.size() == 0)
		return;

	for(iterm = 0; iterm < _nTerms; ++iterm)
		terms_d[iterm] = _terms[iterm];

	interpolate_linear(map, terms_d, coef);
	_eod_val = _eod_saved;
	for(iterm = 0; iterm < _rets_saved.size2(); ++iterm) {
		_eod_val[iterm] = _eod_saved[iterm] * coef[iterm];
	}
}
void Factor::scale_rets(std::map<double,double>& map) {
	unsigned int iterm;
	std::vector<double> terms_d(_nTerms), coef(_nTerms);

	if(_rets_saved.size2() == 0)
		return;

	for(iterm = 0; iterm < _nTerms; ++iterm)
		terms_d[iterm] = _terms[iterm];

	interpolate_linear(map, terms_d, coef);

	_sim_rets = _rets_saved;
	for(iterm = 0; iterm < _rets_saved.size2(); ++iterm) 
		column(_sim_rets, iterm) = column(_rets_saved, iterm) * coef[iterm];
	
}
void Factor::save_eod() {
	if(_eod_val.size() > 0)
		_eod_saved = _eod_val;
}
void Factor::save_rets() {
	if(_sim_rets.size2() > 0)
		_rets_saved = _sim_rets; 
}
void Factor::reset_eod() {
	if(_eod_saved.size() == _nTerms ) 
		_eod_val = _eod_saved;
}
void Factor::reset_rets() {
	if(_nTerms > 0 && _rets_saved.size2() == _nTerms) 
		_sim_rets = _rets_saved; 
}
void Factor::assignEODVal(std::map<double, double>& mp, const FactSet& pos) {
	std::vector<boost::shared_ptr<Factor> >::const_iterator jtfac;
	std::vector<double> arg, val;
	unsigned int iterm;

	if(_mapping_type == Proxied) { // Proxiied factor tenors are identical to proxy tenors
		jtfac = pos.findProxyFact(_proxy_name);
		if(jtfac != pos.getFactors().end()) {
			_terms = (*jtfac)->getTerms();
			_nTerms = (*jtfac)->getNterms();
		}
	}
	arg.resize(_nTerms);
	for(iterm = 0; iterm < _nTerms; ++iterm)
		arg[iterm] = _terms[iterm];
	interpolate_fwd_linear(mp, arg, _eod_val);
	
	mp.clear();
}
void Factor::setRegBasedOnCur(const std::map<std::string, std::string>& cur_2_reg_map) {
	std::map<std::string, std::string>::const_iterator itmap(cur_2_reg_map.find(_curr));
	std::string defaultReg("Western Europe");
	if(itmap != cur_2_reg_map.end())
		_reg = itmap->second;
	else
		_reg = defaultReg;
}
void Factor::apply_CCAR_info(const std::vector<boost::shared_ptr<CCAR_shift> >& ccar_shifts, unsigned int i_proj, const std::map<std::pair<std::string, std::string>, std::string>& mp, 
	const std::map<std::string, std::string>& reg_2_cur_mp, const std::map<std::string, std::string>& cur_2_reg_map) {
	
	//std::vector<boost::shared_ptr<Factor> >::const_iterator itfac;
	std::string core_fact(""), core_cur(""), type_name("");
	//, core_vol_fact(":SP500_C00000117_EQTV");
	
	std::map<std::pair<std::string, std::string>, std::string>::const_iterator itmap;
	std::map<double, double> mpi;
	std::map<std::string, std::string>::const_iterator itmp_reg_2_cur;
	std::pair<std::string, std::string> pair;
	std::vector<boost::shared_ptr<CCAR_shift> >::const_iterator it_shifts;
	std::vector<double> tenors, vals, tt;
	matrix<double> shifts;
	unsigned int nterms, j;

	if(_type == MBS && (_s_name == ":REFI_DAILY" || _s_name == ":TURNOVER_DAILY"))
		return;

	// From core CCAR currency corresponding to the region
	itmp_reg_2_cur = reg_2_cur_mp.find(_reg);
	if(itmp_reg_2_cur == reg_2_cur_mp.end())
		core_cur = "EUR";
	else
		core_cur = itmp_reg_2_cur->second;

	type_name = factTypeNames[_type];

	// Find core CCAR factor name
	pair = std::make_pair<std::string, std::string>(type_name, core_cur);
	itmap = mp.find(pair);
	if(itmap == mp.end())
		return;
	core_fact = itmap->second;
	tt.resize(_nTerms);
	for(j = 0; j < _nTerms; ++j)
		tt[j] = _terms[j];

	// Find core CCAR shifts, change levels 
	for(it_shifts = ccar_shifts.begin(); it_shifts != ccar_shifts.end(); ++it_shifts) {
		if((*it_shifts)->getName() == core_fact ) {
			// Interpolate shifts
			tenors = (*it_shifts)->getTerms();
			nterms = tenors.size();
			shifts = (*it_shifts)->getShifts();
			vals.resize(nterms);
			mpi.clear();
			for(j = 0; j < nterms; ++j)
				mpi.insert(std::make_pair<double, double>(tenors[j], shifts(j,  i_proj)));
				
			interpolate_linear(mpi, tt, vals);
				
			if((*it_shifts)->getIsLevel()) {
				if((*it_shifts)->getIsMult()) 
					for(j = 0; j < _nTerms; ++j)
						_eod_val[j] *= vals[j];
				else
					for(j = 0; j < _nTerms; ++j)
						_eod_val[j] += vals[j];
			}
			break;
		}
	}
	// Change returns
	for(it_shifts = ccar_shifts.begin(); it_shifts != ccar_shifts.end(); ++it_shifts) {
		if(!(*it_shifts)->getIsLevel()) {
			// Interpolate shifts
			tenors = (*it_shifts)->getTerms();
			nterms = tenors.size();
			shifts = (*it_shifts)->getShifts();
			vals.resize(nterms);
			mpi.clear();
			for(j = 0; j < nterms; ++j)
				mpi.insert(std::make_pair<double, double>(tenors[j], shifts(j,  i_proj)));
							
			interpolate_linear(mpi, tt, vals);
				
			if(!(*it_shifts)->getIsMult()) // Additive shift in vol: new vola/ old vol = (vol+vals)/vol
				for(j = 0; j < _nTerms; ++j)
					vals[j] = 1.0 + vals[j] / _vol[j];

			for(j = 0; j < _nTerms; ++j)
				column(_sim_rets, j) *= vals[j];
				
			break;
		}
	}
}

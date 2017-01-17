//#include <boost/date_time/gregorian/gregorian.hpp>
//#include <boost/date_time/gregorian/greg_month.hpp>
//using namespace boost::gregorian;

/*
int alignDates(vector<boost::shared_ptr<Factor> >& factors, vector<boost::shared_ptr<Factor> >& factorsGood, vector<boost::shared_ptr<Factor> >& factorsWithFewDates,vector<date>& dts_common) {
	int maxDates(3600); // about 10 years daily
	int minDatesCount(480); // for regular VaR based on 2 years of data. Should be 240 for Stressed VaR based on 1 year. To generalize Temp Tanya
	vector<date> dts, dts_res;
	int nfactors, i, j, iterm, nterms, ndates, idate;
	vector<boost::shared_ptr<Factor> >::iterator itfac;
	vector<date>::iterator it_res;
	dts_res.resize(maxDates);
	std::ofstream log_file;
	std::map<date, vector<double> >::iterator it;
	vector<double> prc;
	log_file.open("G:\\VAR attribution\\20140224\\Output\\log_histVar.csv");

	nfactors = factors.size();
	i = 0;
	itfac = factors.begin();
	while( (*itfac)->getHrDates().size() < minDatesCount) {
		log_file << "alignDates: " << (*itfac)->getScenName() << ": too few dates " << "," << dts_res.size() << endl;
		factorsWithFewDates.push_back(std::move((*itfac)));
		++itfac;
		++i;
	}

	// Sort factors that have sufficiently long history and positive prices from the rest
	if(itfac != factors.end()) {
		dts_common = (*itfac)->getHrDates();
		++itfac;
		++i;
		while(itfac != factors.end()) {
		
			dts = (*itfac)->getHrDates();
			dts_res.resize(maxDates);
			it_res = set_intersection(dts_common.begin(), dts_common.end(), dts.begin(), dts.end(), dts_res.begin());
			dts_res.resize(it_res - dts_res.begin());
			if(dts_res.size() < minDatesCount) { // exclude this factor 
				log_file << "alignDates: " << (*itfac)->getScenName() << ": too few dates " << "," << dts_res.size() << endl;
				factorsWithFewDates.push_back(std::move((*itfac)));
			}
			else	{
				factorsGood.push_back(std::move((*itfac)));
				dts_common.assign(dts_res.begin(), dts_res.end());
				j = dts_common.size();
			}
			++i;
			itfac = factors.begin() + i;
		}
	}
	ndates = dts_common.size();
	nfactors = factors.size();
	
	for(itfac=factorsGood.begin(); itfac != factorsGood.end(); ++itfac) {
		nterms = (*itfac)->getNterms();
		
		//matrix<double>* matr = new matrix<double>(ndates, nterms);
		//matrix<double>* rets = new matrix<double>((ndates-1), nterms);
		(*itfac)->_histData.resize(ndates, nterms);
		(*itfac)->_rets.resize(ndates - 1, nterms);
		
		for(idate = 0; idate < ndates; ++idate) {
			it = (*itfac)->_hrDataMap.find(dts_common[idate]);

			if(it!=(*itfac)->_hrDataMap.end() ) {
				prc = (*it).second;
				if(prc.size() == nterms) {

					for(iterm = 0; iterm < nterms; ++iterm)
						(*itfac)->_histData(idate, iterm) = prc[iterm]; 
						//(*matr)(idate, iterm) = prc[iterm];
					
					if(idate > 0)
						for(iterm = 0; iterm < nterms; ++iterm) {
							if(abs((*itfac)->_histData(idate - 1, iterm)) > 1E-10)
								(*itfac)->_rets(idate-1, iterm) = log(prc[iterm] / (*itfac)->_histData(idate - 1, iterm)); 
							else
								log_file << "alignDates: zero historical price " << (*itfac)->getScenName() << endl;
						}
					
					//(*rets)(idate-1, iterm) = log(prc[iterm] / (*matr)(idate - 1, iterm));
					if(idate == ndates -1)
						(*itfac)->setSpot(prc);
				}
				else
					log_file << "alignDates: number of terms mismatch for " << (*itfac)->getScenName() << endl;
			}
		}
		
		//(*itfac)->setHistData(matr);
		//(*itfac)->setRets(rets);
	}
	
	log_file.close();
	return 0;
}
*/
/*
int read_FX_maps(vector<boost::shared_ptr<Factor> >& s_facts, std::map<string, vector<string> >& reg_2_cur_map) {
	// For each region there are a number of currencies that can be used in regression for the factors with insufficient data
	// For some of this currencies there may not be positions. In this case create factor with zero sensitivity
	// Also a map to hold the mapping is created here
	string buf(""), factName, cur, subType, reg, subtype("");
	vector<string> values;
	vector<double> sens(1, 0);
	int errCode(-1002), importance, nterms(0);
	vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::map<string, vector<string> >::iterator itmap;

	string fx_reg_2_cur_fn("G:\\VAR attribution\\20140224\\Maps\\FX_map_reg_2_cur.csv");
	std::ifstream f_reg_2_cur;
	f_reg_2_cur.open(fx_reg_2_cur_fn);
	
	using boost::is_any_of; 
	getline(f_reg_2_cur, buf);
	boost::algorithm::split(values, buf, is_any_of(",")); 
	if(values.size() != 3)
		return errCode;

	while(getline(f_reg_2_cur, buf)) {
		boost::algorithm::split(values, buf, is_any_of(",")); 
		if(values.size() != 3)
			return errCode;
		reg = values[0].c_str();
		cur = values[1].c_str();
		importance = atoi(values[2].c_str());

		itfac = findFactor(FX_Spot, subType, cur, s_facts);
		if(itfac != s_facts.end()) {
			(*itfac)->setImportance(importance);

		}	
		else { // If therse no such factor yet
			boost::shared_ptr<Factor> ptr(new Factor);
			vector<int> terms(1,0);
			ptr->setAttributes(FX_Spot, subType, cur, reg, sens, terms, importance);
			s_facts.push_back(std::move(ptr));
		}
		try {
			reg_2_cur_map.at(reg).push_back(cur);
		}
		catch(std::out_of_range) {
			vector<string> tmp;
			tmp.push_back(cur);
			reg_2_cur_map.insert(std::pair<string, vector<string> >(reg, tmp));
		}
	}

	f_reg_2_cur.close();

	return 0;
}
*/
/*
void gramm_shmidt_orth(matrix<double>& mat) { 
// Orthogonalize rows of mat
	matrix<double> mcp(mat);
	double ss;
	std::vector<double> wts, rw;
	int i, j, k, n(mat.size1()), m(mat.size2());

	for(j = 0; j < n; ++j) {
		rw = row(mat, j);
		
		wts.resize(j + 1);
		for(k = 0; k < j; ++k) {
			wts[k] = inner_product(row(mat, j).begin(), row(mat, j).end(), row(mat, k).begin(), 0.0);


		for(k = 0; k < j; ++k)	
			rw -= wts[k] * row(mat, k);

		ss = inner_product(rw.begin(), rw.end(), rw.begin(), 0.0);
		ss = 1.0/sqrt(ss);
		for(i = 0; i < m; ++i) 
			mat(j, i) = ss * rw[i];
	}
}
*/
/*
void gramm_shmidt_orthog(matrix<double>& mat) { 
// Orthogonalize rows of mat
	matrix<double> mcp(mat);
	double ss;
	std::vector<double> wts, rw1, rw2;
	int i, j, k, n(mat.size1()), m(mat.size2());

	for(j = 0; j < n; ++j) {
		for(i = 0; i < m; ++i)
			mat(j, i) = mcp(j, i);

		wts.resize(j + 1);
		for(k = 0; k < j; ++k) 
			wts[k] = rows_scaler_prod(mat, j, k);

		for(k = 0; k < j; ++k)	{

			//for(i = 0; i < m; ++i)
			//	mat(j, i) -= wts[k] * mat(k, i);
		}
		ss = rows_scaler_prod(mat, j, j);
		ss = 1.0/sqrt(ss);
		for(i = 0; i < m; ++i) 
			mat(j, i) *= ss;
	}
}
*/
/*
double rows_scaler_prod(matrix<double>& mat, int n1, int n2) {
	double ss(0.0);
	int i, m(mat.size2());

	for(i = 0; i < m; ++i)
		ss += mat(n1, i) * mat(n2, i);

	ss /= m;

	return ss;
}
*/
	// Temp Tanya: for now the file with equity spot names is not available. For testing purpose only make up spot rf names from equity vol risk factor names.
	// For example: 
	//:AAPL_C03783310_EQTV		C03783310_NASDAQ_AAPLSpot
	// :AA_C01381710_EQTV		C01381710_NYSE_AASpot,
	// :ABBV_C00287Y10_EQTV		C00287Y10_NYSE_ABBVSpot
	//:ABT_C00282410_EQTV		C00282410_NYSE_ABTSpot
	// :ABX_C06790110_EQTV		C06790110_NYSE_ABXSpot
	/*
	string eq_spot_names[] = {"C03783310_NASDAQ_AAPLSpot","C01381710_NYSE_AASpot", "C00287Y10_NYSE_ABBVSpot", "C00282410_NYSE_ABTSpot", "C06790110_NYSE_ABXSpot"};
	std::vector<string> eqt_sp(eq_spot_names, eq_spot_names + sizeof(eq_spot_names)/sizeof(eq_spot_names[0]));
	unsigned int eqt_n(eqt_sp.size());
	for(i = 0; i < eqt_n; ++i)
		names.push_back(eqt_sp[i]);
	*/
	
	/*
int read_FX_maps(string map_fn, std::map<string, std::vector<string> >& reg_2_cur_map) {
	// For each region there are a number of currencies that can be used in regression for the factors with insufficient data
	// Such map is created here
	string buf(""), factName, cur, subType, reg, subtype("");
	std::vector<string> values;
	std::vector<double> sens(1, 0);
	int errCode(-1002), nterms(0);
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::map<string, std::vector<string> >::iterator itmap;

	std::ifstream f_reg_2_cur;
	f_reg_2_cur.open(map_fn);
	
	using boost::is_any_of; 
	getline(f_reg_2_cur, buf);
	boost::algorithm::split(values, buf, is_any_of(",")); 
	if(values.size() != 3)
		return errCode;

	while(getline(f_reg_2_cur, buf)) {
		boost::algorithm::split(values, buf, is_any_of(",")); 
		if(values.size() != 3)
			return errCode;
		reg = values[0].c_str();
		cur = values[1].c_str();
		// importance = atoi(values[2].c_str()); // not used yet

		try {
			reg_2_cur_map.at(reg).push_back(cur);
		}
		catch(std::out_of_range) {
			std::vector<string> tmp;
			tmp.push_back(cur);
			reg_2_cur_map.insert(std::pair<string, std::vector<string> >(reg, tmp));
		}
	}

	f_reg_2_cur.close();

	return 0;
}
*/
/*
bool Factor::mapIt(double defaultIdioVol, unsigned int minNumDatesForRegr, std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Spot, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Spot, 
	std::vector<boost::shared_ptr<Factor> >& f_IR_Curve, 	std::vector<boost::shared_ptr<Factor> >& f_FX_Spot, std::vector<boost::shared_ptr<Factor> >&f_FX_Forward, std::vector<date>& dts, 
	std::map<std::string, std::vector<std::string> >& reg_2_cur_map) {

	bool success(false);
	unsigned int j, i, ndates, nrets, nfacs, ifac, iterm, theTermIdx, theTerm;

	std::vector<boost::shared_ptr<Factor> >:: iterator jt;
	std::vector<boost::shared_ptr<Factor> > factors_4_regression;
	std::vector<std::string> currencies_4_mapping;
	std::vector<std::string>::iterator itcur;
	std::vector<date> hrDates;
	std::vector<date>::iterator itdt, jtdt;
	std::vector<unsigned int> terms;
	std::vector<double> prc, ave, R2, thisFactorData, r2;
	std::string subType("");
	double idio_vol(defaultIdioVol);
	std::map<date, std::vector<double> >::iterator itmap;
	std::map<date, std::vector<double> > hrDataMap;
	matrix<double> rets, rets_all, histData, regr_coef;
		
	currencies_4_mapping = reg_2_cur_map.at(_reg);
	itdt = set_intersection(_hrDates.begin(), _hrDates.end(), dts.begin(), dts.end(), _hrDates.begin());
	_hrDates.resize(itdt - _hrDates.begin());

	if(_h_name == "C49926D10_NYSE_KNSpot" )
		int tt = 0;

	if(_hrDates.size() < minNumDatesForRegr || _nTerms == 0) {
		setRetsVols(defaultIdioVol, minNumDatesForRegr, _hrDates);
		return false;
	}
	ndates = _hrDates.size();
	nrets = ndates - 1;

	// Find set of factors to regress against
	switch(_type) {
		case IR_Curve:	// Regress against data for the same region
			for(itcur = currencies_4_mapping.begin(); itcur != currencies_4_mapping.end(); ++itcur) {
				jt = findFactor(_type, _subType, *itcur, f_IR_Curve);
				if(jt != f_IR_Curve.end())
					factors_4_regression.push_back((*jt));
			}
		case FX_Forward:
			for(itcur = currencies_4_mapping.begin(); itcur != currencies_4_mapping.end(); ++itcur) {
				jt = findFactor(_type, _subType, *itcur, f_FX_Forward);
				if(jt != f_FX_Forward.end())
					factors_4_regression.push_back((*jt));
			}
			break;

		case FX_Spot:
			factors_4_regression = f_PCA_FX_Spot;
			
			for(itcur = currencies_4_mapping.begin(); itcur != currencies_4_mapping.end(); ++itcur) {
				jt = findFactor(_type, _subType, *itcur, f_FX_Spot);
				if(jt != f_FX_Spot.end())
					factors_4_regression.push_back((*jt));
			}
			
			break;
		case Eqt_Spot: // Regress agains indices
			factors_4_regression = f_PCA_Eqt_Spot;
			break;

		case CDS:
			idio_vol = defaultIdioVol;
			break;
		case FX_Pair_Vol:
			idio_vol = defaultIdioVol;
			break;
		case EQ_Vol: 
			idio_vol = defaultIdioVol;
			success = true;
			break;
		case IRCap_Vol:
			idio_vol = defaultIdioVol;
			success = true;
			break;
		case IRSwaption_Vol:
			idio_vol = defaultIdioVol;
			success = true;
			break;
		default:
			idio_vol = defaultIdioVol;
			break;
		// Turnover, Refi, MBS, factor_GB, 
	}
	nfacs = factors_4_regression.size();
	if(nfacs == 0) {
		setRetsVols(defaultIdioVol, minNumDatesForRegr, _hrDates);
		return false;
	}
	
	// Collect historical data for regression
	setRetsVols(defaultIdioVol, minNumDatesForRegr, dts);
	nrets = _hrDates.size() - 1;
	if(nrets < minNumDatesForRegr)
		return false;

	// If factor is mapped the number of terms is set to 1.
	theTermIdx = _nTerms * 0.5;
	theTerm = _terms[theTermIdx];
	_theTermForMapping = theTermIdx;

	matrix<double> rets_this(nrets, 1);
	_R2.resize(1);
	_regr_coef.clear();
	_regr_coef.resize(1, nfacs);

	column(rets_this, 0) = column(_rets, theTermIdx);
	thisFactorData.resize(nrets + 1);
	rets_all.resize(nrets, nfacs);

	ifac = 0;
	for(jt = factors_4_regression.begin(); jt != factors_4_regression.end(); ++jt) {
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
		hrDataMap = (*jt)->_hrDataMap;

		for(j = 0; j < _hrDates.size(); ++j) {
			prc = hrDataMap.find(_hrDates[j])->second;
			thisFactorData[j] = prc[jj];
			if(j > 0)
				rets_all(j - 1, ifac) = log(thisFactorData[j] / thisFactorData[j - 1]);
		}
		_mapped_fact.push_back(*jt);
		_mapped_terms.push_back(jj);
		++ifac;
	}

	// Regression
	success = regress(rets_all, rets_this, regr_coef, _resid_vol, r2);
	if(success) {
		iterm = 0;
		_R2[iterm] = r2[0];
		for(ifac = 0; ifac < nfacs; ++ifac)
			_regr_coef(iterm, ifac) = regr_coef(ifac,0);
	}
	
	return success;
}

double volatility(std::vector<double>& vec) {
	double vol(0.0);
	unsigned int n(vec.size()), n1(n-1), i;
	std::vector<double> rets(n1);
	if(n > 2) {
		for(i = 1; i < n; ++i) 
			rets[i-1] = log(vec[i] / vec[i-1]);
		vol = stdev(rets);
	}
	return vol;
}
double average(std::vector<double>& vec) {
	double ave(0.0);
	unsigned int n(vec.size()), i;
	if(n > 0) {
		for(i = 0; i < n; ++i) 
			ave += vec[i];

		ave /= n;
	}
	return ave;
}

*/
/* Setup random numbers generator
	boost::mt19937                     gener(1);
    boost::normal_distribution<double> normal_dist;   // Normal Distribution
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal_dist);
	rng.engine().seed(1);
	int nsim(1000), nfact(4);
	matrix<double> randoms(nsim, 4);
	double vv(0);

	for(int k = 0; k < nsim; ++k) 
			for(i = 0; i < nfact; ++i) 
				randoms(k, i) = rng();
	if(1) {  
		gramm_shmidt_orth(randoms); // Orthogonalize columns
		vv = orthog_check(randoms);
	}
	*/
	/* Test
	matrix<double> X(10, 2), Y(10, 3), coef;
	std::vector<double> vols, R2;
	int i, j, nn(10), mm(2), kk(3);
	for(i = 0; i < nn; ++i) {
		for(j = 0; j < mm; ++j) 
			X(i, j) = i + (1+ 0.01 *(i+j));
		for(j = 0; j < kk; ++j)
			Y(i, j) = i +0.0001* j;
	}
	regress(X, Y, coef, vols, R2);
	*/
/*
int genShiftScen(std::vector<boost::shared_ptr<Factor> >& factors, string shifts_dir, string asOfDate) {
	
	// Set of market risk factors will be based on the latest historical data;
	// Shifts will be written to a file plus info in  fact_info file
	std::vector<unsigned int> iterms, jterms;
	std::ofstream out_file;
	const string line1("Scenario Set,Scenario Set Name,Scenario Name,Scenario Probability,Scenario Color,Scenario Variable,Scenario Start Time,Scenario Attribute,Time Evolution,Time Evolution to Trigger,Trigger Holder,	Scenario Shift Rule,Scenario Type,Scenario Replacement Value,");
	unsigned int i, j, m, itp, k(1);

	// Generate shifts by factor type
	double shiftVal, shft1bp(0.0001);

	const string scen_color("green"), scenSetName("msmc 2000"), scen_name("Mmc2000"), ScenarioStartTime(asOfDate), TimeEvolution("Constant"), TimeEvolutionToTrigger("Constant"),
		TriggerHolder("@Standard()"), ScenShiftRule("Trigger Time"), ScenType("non-parallel shift"), ScenReplVal("Term"), Attr14(asOfDate);

	std::vector<boost::shared_ptr<Factor> >::iterator it, jt;
	FactType iType;

	for(itp = 0; itp < Factor::factTypeNames.size(); ++itp) {
		iType = FactType(itp);

		out_file.open(shifts_dir + Factor::factTypeNames[itp] + ".csv");
		out_file << line1 << endl;
		k = 1;
		for(it = factors.begin(); it != factors.end(); ++it) {
				if((*it)->getType() == iType ) { 
					iterms = (*it)->getTerms();
					if(iterms.size() > 0 && (*it)->getScenName() != "") {
						for(i = 0; i < iterms.size(); ++i) { // new scenario
							// In each scenario single factor is perturbed
							if(k == 1)
								out_file << "," << scenSetName << "," << scen_name << k << "," << 0.001;
							else
								out_file << ",," << scen_name << k << "," << 0.001;
							m = 0;
							for(jt = factors.begin(); jt != factors.end(); ++jt) {
								if((*jt)->getType() == iType && (*jt)->getTerms().size() > 0 && (*jt)->getScenName() != "") { // new factor
									if(m != 0)
										out_file << ",,,";
									out_file << "," << scen_color << "," << (*jt)->getScenName() << "," << ScenarioStartTime << ",," << TimeEvolution << "," << TimeEvolutionToTrigger << "," << 
										TriggerHolder << "," << ScenShiftRule << "," << ScenType << "," << ScenReplVal << "," << Attr14 << endl;
								
									jterms = (*jt)->getTerms();
									for(j = 0; j < jterms.size(); ++j) {
										shiftVal = 0.0;
										if(k != 0 && it == jt && i == j)
											shiftVal = shft1bp;
						
										out_file << ",,,,,,,,,,,,," << jterms[j] << "," << shiftVal << endl;
									}
								++m;
								} // endif
							} // end for
							++k;
						} // end for
					} // end if(terms.size() > 0)
				}
			}
			// Add base scenario
			m = 0;
			shiftVal = 0.0;
			out_file << ",," << scen_name << k << "," << 0.001;
			for(jt = factors.begin(); jt != factors.end(); ++jt) {
				if((*jt)->getType() == iType) { // new factor
					out_file << ",";
					if(m != 0)
						out_file << ",,,";
					out_file << scen_color << "," << (*jt)->getScenName() << "," << ScenarioStartTime << ",," << TimeEvolution << "," << TimeEvolutionToTrigger << "," << 
						TriggerHolder << "," << ScenShiftRule << "," << ScenType << "," << ScenReplVal << "," << Attr14 << endl;
					jterms = (*jt)->getTerms();
					for(j = 0; j < jterms.size(); ++j) 
						out_file << ",,,,,,,,,,,,," << jterms[j] << "," << shiftVal << endl;
					++m;
			} // endif
		} // end for
		out_file.close();
	}
 	return 0;
}
int setSensitivities(string dir, std::vector<boost::shared_ptr<Factor> >& factors, std::map<string, string>& cur_2_reg_map) {
	// Read sensitivities from files, they correspond to risk factors in the same order as in shift scenarios
	FactType iType;
	string buf(""), tmp, cur, reg;
	unsigned int i, itp, scenNum(0), iscen(0);
	string line1;
	std::ifstream sens_file;
	std::vector<unsigned int> iterms;
	std::vector<double> vals;
	std::vector<string> values;
	std::vector<boost::shared_ptr<Factor> >::iterator it;
	std::map<string, string>::iterator itmap;

	for(itp = 0; itp < Factor::factTypeNames.size(); ++itp) {
	  iType = FactType(itp);
	  if(iType == FX_Spot || iType == FX_Forward || iType == IR_Curve || iType == FX_Pair_Vol) { // temp Tanya
		//string shift_fn = shifts_dir + "shift_" + Factor::factTypeNames[itp] + ".csv";
		string shift_fn = dir + "Shifts\\" + Factor::factTypeNames[itp] + ".csv";
		string sens_fn = dir + "Sens\\" + Factor::factTypeNames[itp] + ".csv";
		std::vector<boost::shared_ptr<Factor> > factors_this_type;
		
		readScenarioFile(shift_fn, factors_this_type, line1);
		sens_file.open(sens_fn);

		if(sens_file.is_open()) {
			getline(sens_file, buf);
			scenNum = atoi(buf.c_str());
			iscen = 0;
			for(it = factors_this_type.begin(); it != factors_this_type.end(); ++it) {
				iterms = (*it)->getTerms();
				for(i = 0; i < iterms.size(); ++i)  {
					if(getline(sens_file, buf)) {
						boost::algorithm::split(values, buf, is_any_of(" "));
						vals.push_back(atof(values[0].c_str()));
						
					}
				}
				++iscen;
				(*it)->setSens(vals);
				vals.clear();
				factors.push_back(std::move(*it));
			}
			sens_file.close();
			
		}
	  }
	} 
	return 0;
}
*/
/*
// Output vars by factor
		ofstream var_out;
		var_out.open(dir_var + "var_report_all.csv");
		var_out << "All factors VAR:," << var_MC << endl;
		var_out << "Term, Name, Mapping Type, Asset Class, Currency, region, EOD value, Hist Vol, Ratio of Idiosyncr and Sstematic variances, Delta, Gamma, VaR by term, VaR" << endl;
		for(i = 0; i != h_factors.size(); ++i)
			h_factors[i]->print_vars(var_out);
		var_out.close();

		var_MC = calc_MC_Var(f_noGaps, percentile);
		var_out.open(dir_var + "var_report_noGaps.csv");
		var_out << "No Gaps factors VAR:," << var_MC << endl;
		var_out << "Term, Name, Mapping Type, Asset Class, Currency, region, EOD value, Hist Vol, Ratio of Idiosyncr and Sstematic variances, Delta, Gamma, VaR by term, VaR" << endl;
		for(i = 0; i != f_noGaps.size(); ++i)
			f_noGaps[i]->print_vars(var_out);
		var_out.close();

		var_MC = calc_MC_Var(f_Idiosyncr, percentile);
		var_out.open(dir_var + "var_report_Idio.csv");
		var_out << "Idio factors VAR:," << var_MC << endl;
		var_out << "Term, Name, Mapping Type, Asset Class, Currency, region, EOD value, Hist Vol, Ratio of Idiosyncr and Sstematic variances, Delta, Gamma, VaR by term, VaR" << endl;
		for(i = 0; i != f_Idiosyncr.size(); ++i)
			f_Idiosyncr[i]->print_vars(var_out);
		var_out.close();

		var_MC = calc_MC_Var(f_Mapped, percentile);
		var_out.open(dir_var + "var_report_Mapped.csv");
		var_out << "Mapped factors VAR:," << var_MC << endl;
		var_out << "Term, Name, Mapping Type, Asset Class, Currency, region, EOD value, Hist Vol, Ratio of Idiosyncr and Sstematic variances, Delta, Gamma, VaR by term, VaR" << endl;
		for(i = 0; i != f_Mapped.size(); ++i)
			f_Mapped[i]->print_vars(var_out);
		var_out.close();
*/
/*
int _stdcall process(LPCSTR stockmap_fn, LPCSTR varscen_fn, LPCSTR factorCoef_fn, LPCSTR RAT_fn, LPCSTR volFactList_fn, long nScen_long, LPCSTR date, LPCSTR vars_fn, LPCSTR pfls_fn)
{
	// Reads the following files related to Equities: stockmap, varscen, factorCoef, RAT, volFactList
	// Calculates VaR per position based on varscen and VaRs per portfolios
	const int nPC(18);
	unsigned int factorCoefMax(50000);
	// In stockmap.dat the equity price name for S&P is C00000117_Index_SPALN, the volatility name is SP500_C00000117
	// In RAT the ticker name SP500
	const string volNameSPALN("SP500");
	std::vector<string> pcNames, volFacNames;
	unsigned int i, j, k, icount(0), nScen(nScen_long);
	string ss, ff(":factorGB_");
	for(i = 0; i < nPC; ++i) {
		string nn = std::to_string(static_cast<long long>(i));
		ss = ff + nn;
		pcNames.push_back(ss);
	}

	using namespace boost::numeric::ublas;
	using boost::is_any_of; 
	
	std::vector<string> values, sub_values;
	string buf("");
	string curName(""), curName2("");
	
	std::vector<string> eqNames;
	matrix<double> eqScen;
	read_scen(varscen_fn, nScen, eqNames, eqScen);

	// Read volatility factor names
	std::ifstream volfacnames_file;
	volfacnames_file.open(volFactList_fn);
	icount = 0;
	while(getline(volfacnames_file, buf)) {
		boost::algorithm::split(values, buf, is_any_of(","));
		volFacNames.push_back(values[0]);
	}

	// Separate volatilities simulation from 18 factors simulations
	matrix<double> eqFacScen(nScen, nPC);
	matrix<double> eqVolScen(nScen, volFacNames.size());
	icount = 0;
	for(i = 0; i < eqNames.size(); ++i) {
		j = 0;
		while(j < nPC && eqNames[i].compare(pcNames[j]))
			j++;
		if(j < nPC)
			for(k=0; k < nScen; ++k)
				eqFacScen(k, j) = eqScen(k, i);
		else
			j = 0;
			while(j < volFacNames.size() && eqNames[i].compare(volFacNames[j]))
				++j;
			if(j < volFacNames.size()) {
				for(k=0; k < nScen; ++k)
					eqVolScen(k, icount) = eqScen(k, i);
				++icount;
			}
	}
	
// read stock map file
	std::ifstream stockmap_file;
	std::map<string, int> stockMap;
	std::map<string, string> stockCurMap;
	std::map<string, int> stockNumTenorsMap;
	stockmap_file.open(stockmap_fn);
	std::ofstream nonStandardNamesStockmapDat;
	nonStandardNamesStockmapDat.open("G:\\Var Attribution\\temp\\nonStandardNameStockmap.csv");
	icount = 0;
	int ct(0), ind(0);
	while(getline(stockmap_file, buf)) {
		// File structure: "ACS_C00819010	0	USD	1" - equityName_cusip interger currency tenor
		boost::algorithm::split(values, buf, is_any_of("|"));	
		
		if(curName == values[0])
			++ct;
		else {
			boost::algorithm::split(sub_values, values[0], is_any_of("_"));
			ind = atoi(values[1].c_str());
			if(sub_values.size()==3) {
				curName2 = sub_values[0] + "_" + sub_values[2];
				stockMap.insert(std::pair<string, int >(curName2, ind));
			}
			else
				nonStandardNamesStockmapDat << values[0] << "," << ind << "," << sub_values.size() << endl;

			curName = values[0];
			ct = 0;
		}
		++icount;
	}
	nonStandardNamesStockmapDat.close();

	// read equities factor loadings
	matrix<double> factors(factorCoefMax, nPC);
	std::vector<double> idioRisk;
	std::ifstream file;
	file.open(factorCoef_fn);

	icount=0;
	j=0;
	while(getline(file, buf)) {
		boost::algorithm::split(values, buf, is_any_of(","));
		if(icount > 1) {
			idioRisk.push_back(atof(values[10].c_str()));
			for(int i = 0; i < nPC; ++i) {
				factors(icount - 2, i) = atof(values[i+11].c_str());
				
				if(abs(factors(icount - 2, i))>1000) {
					double v(fabs(factors(icount - 2, i)));
					v = v;
				}
			}
		}
		++icount;
	}

	// read RAT
	string sedol, cusip, secId, eqName, pfl, cur, maturity, discCrv, instrType, ticker, exchange;
	double spot, quantity, not, theoVal, delta, gamma, BpsUp, BpsDn, vegaUp, vegaDn, beta, var, calcVar, idio; 
	
	std::ifstream rat_file;
	//string rat_file_path("G:\\Var Attribution\\SlnInput\\RAT.csv");
	//string rat_file_path("G:\\Var Attribution\\SlnInput\\RAT_MUREX_20130319.txt");
	rat_file.open(RAT_fn);
	
	std::vector<EqPosPtr> pos;
	std::vector<portfolio> pfls;
	std::vector<string> tickerCusip;
	int notFoundVol(0), notFoundEq(0);
	std::vector<string> portfolios;

	std::ofstream file_out;
	file_out.open(vars_fn);
	file_out << "secName" << "," << "portfoli" << "," << "spot" << "," << "not" << "," << "delta" << "," << "gamma" << "," << "vegaDn" << "," << "vegaUp" << "," << "eqInd" << "," << 
		"eqVolInd" << "," << "var" << "," << "calcVar" << "," << "idio";
	file_out << endl;
	while(getline(rat_file, buf)) {
		if(buf.find("SourceID")!=string::npos) 
			continue;
                                                                                                                                 
		boost::algorithm::split(values, buf, is_any_of(","));
		sedol = values[51];
		cusip = values[50];
		secId = values[49];
		ticker = values[53];
		
		pfl = values[52];
		portfolios.push_back(pfl);
		cur = values[25];
		maturity = values[15];
		discCrv = values[30];
		instrType = values[9];
		spot = atof(values[37].c_str());
		quantity = atof(values[20].c_str());
		not = atof(values[19].c_str());
		theoVal = atof(values[22].c_str());
		delta = atof(values[24].c_str());
		gamma = atof(values[57].c_str());
		BpsUp = atof(values[31].c_str());
		BpsDn = atof(values[32].c_str());
		vegaUp = atof(values[33].c_str());
		vegaDn = atof(values[34].c_str());
		beta = atof(values[55].c_str());
		var = atof(values[45].c_str());
		exchange = values[8];
		
		string volName, secName;
		if(!ticker.compare("SP500")) {
			volName=":SP500_C" + cusip +"_EQTV";
			secName = "C" + cusip + "_" + "SPALN";
			//secName = "C" + cusip + "_" + exchange + "_" + "SPALN";
		}
		else {
			if(secId == "S") {
				volName=":" + ticker+"_S"+sedol +"_EQTV";
				secName = "S" + sedol + "_" + ticker;
			}
			else {
				volName=":" + ticker+"_C"+cusip + "_EQTV";
				secName = "C" + cusip + "_" + ticker;
			}
		}

		// Find corresponding equity vol index
		std::vector<string>::const_iterator it;
		int eqVolInd(-1);
		j = 0;
		while(j < volFacNames.size() && volName.compare(volFacNames[j]) )
			++j;
		if(j < volFacNames.size()) 
			eqVolInd = j;
		
		// Find corresponding equity index
		std::map<string, int>::const_iterator jt;
		int eqInd(-1);
		
		for(jt = stockMap.begin(); jt!=stockMap.end(); ++jt) {
			if(jt->first.find(secName) != string::npos) {
				eqInd=jt->second;
				break;
			}
		}
		
		if(eqInd <0)
			++notFoundEq;
		if(eqVolInd <0)
			++notFoundVol;
		
		EqPosPtr p(new EqPosition(sedol, cusip, ticker, secId, secName, pfl, eqInd, eqVolInd, cur, maturity, discCrv, instrType, spot, quantity, not, theoVal, delta, gamma, 
			BpsUp, BpsDn, vegaUp, vegaDn, beta, var, factors));

		calcVar = invalidVar;
		idio = 0;
		if(eqInd > 0) {
			if(delta > 0)
				idio = delta *(exp(idioRisk[eqInd] * varQuantile)-1);
			else
				idio = delta *(exp(- idioRisk[eqInd] * varQuantile)-1);
		}
		calcVar = p->calcVar(eqScen, eqFacScen, idioRisk, percentile); 
		p->setCalcVar(calcVar);
		p->setIdioVar(idio);
		file_out << secName << "," << pfl << "," << spot << "," << not << "," << delta  << "," << gamma << "," << vegaDn << "," << vegaUp << "," << eqInd << "," << eqVolInd << "," <<
			var << "," << calcVar << "," << idio;
		file_out << endl;
		pos.push_back(p);

		//vector<portfolio>::const_iterator ip(pfls.begin());
		//while(ip != pfls.end() && pfl.compare(ip->getName()))
		//	++ip;
		
	}
	file_out.close();
	factors.clear();

	// VAR by portfolio
	std::sort(portfolios.begin(), portfolios.end());
	std::vector<string>::iterator ut;
	ut = std::unique(portfolios.begin(), portfolios.end(), equalStr);
	portfolios.resize(std::distance(portfolios.begin(), ut));
	std::ofstream pflReport;
	pflReport.open(pfls_fn);
	pflReport << "portfolio" << "," << "Positions #" << "," << "Delta" << "," << "Gamma" << "," << "Vega" << "," << "deltaGammaVar" << "," << 
		"vegaVar" << "," << "totalVar" << "," << "Idio";
	pflReport << endl;
	for(unsigned int i = 0; i <portfolios.size(); ++i) {
		string pfl(portfolios[i]);
		std::vector<double> res;
		portVarCalc(pfl, pos, eqVolScen, eqFacScen, res);

		pflReport << pfl << ",";
		for(unsigned int j =0; j < res.size(); ++j)
			pflReport << res[j] << ",";
		pflReport << endl;
	}
	pflReport.close();
	return 0;
}
double getVarFromPnl(std::vector<double>& pnl, double perc) {
	int nscen(pnl.size()), intForVar(perc * nscen);
	double var(0);
	if(intForVar > 0) {
		std::nth_element(pnl.begin(), pnl.begin()+intForVar, pnl.end());
		var = *(pnl.begin()+intForVar);
	}
	return var;
}
bool equalStr(string s1, string s2) {
	return(s1==s2);
}
void portVarCalc(string pfl, const std::vector<EqPosPtr>& pos, matrix<double>& eqVolScen, matrix<double>& eqFacScen, 
	std::vector<double>& res) {
	std::vector<EqPosPtr>::const_iterator it;
	double d, eqVar(0), vegaVar(0), totalVar(0), idio(0), delta(0), gamma(0), vega(0), posNum(0) ; 
	std::vector<double> eqPnl, vegaPnl, totalEqPnl, totalVegaPnl, totalPnl;
	unsigned int i, nScen(eqVolScen.size1());

	for(it = pos.begin(); it!= pos.end(); ++it) {
		EqPosPtr p(*it);
		totalPnl.resize(nScen);
		totalEqPnl.resize(nScen);
		totalVegaPnl.resize(nScen);
		if(p->getPfl() == pfl) {
			eqPnl.resize(0);
			vegaPnl.resize(0);
			p->calcEqPnl(eqFacScen, eqPnl);
			p->calcVegaPnl(eqVolScen, vegaPnl);
			d = p->getIdio();
			idio += (d * d);
			delta += p->getDelta();
			gamma += p->getGamma();
			vega += (p->getVegaUp() - p->getVegaDn()) *0.5;
			++posNum;
			if(eqPnl.size() == nScen) 
				for(i = 0; i < eqPnl.size(); ++i)
					totalEqPnl[i] += eqPnl[i] ;
			if(vegaPnl.size() == nScen)
				for(i = 0; i < eqPnl.size(); ++i)
					totalVegaPnl[i] += vegaPnl[i];
		}
	}
	for(i = 0; i < nScen; ++i)
		totalPnl[i] = totalVegaPnl[i] +totalEqPnl[i];
	res.push_back(posNum);
	res.push_back(delta);
	res.push_back(gamma);
	res.push_back(vega);
	eqVar = getVarFromPnl(totalEqPnl, percentile);
	res.push_back(eqVar);
	vegaVar = getVarFromPnl(totalVegaPnl, percentile);
	res.push_back(vegaVar);
	totalVar = getVarFromPnl(totalPnl, percentile);
	res.push_back(totalVar);
	res.push_back(sqrt(idio));
}
*/
/*
int _stdcall parse_scen(LPCSTR inp_fn, LPCSTR out_fn, bool readEqOnly) {
// TT: currently only reading the scenarios for factors that either don't have term structure of have tenor 184 (for spot and vols)
// should be changed
	std::vector<string> scenTerm1;
	std::vector<string> curveName;
	string scenRecord1;
	
	std::ifstream inp_file;
	std::ofstream out_file;
	inp_file.open(inp_fn);
	
	string buf("");
	std::vector<string> values;
	std::vector<string> nameParts;
	string curName, curNameTerm;
	//int unIdentified(0), ic1(0), ic2(0), ic3(0), ic4(0), ic6(0);

	//long count = 0;
	int scenNum(-1), iFactor(0), nFactors(0);
	using namespace boost::numeric::ublas;
	int scenNumMax(1000);
	int numFacsMax(5000);
    matrix<double> scen_mat (scenNumMax, numFacsMax);
	
	while(getline(inp_file, buf)) {
		if (buf.length() == 0)
			continue;
		if(buf.find("base")!=string::npos) 
				break;
		if (buf.find("Scenario Set") == string::npos){ // start new scenario
			if (buf.find("Mmc2000")!=string::npos) {
				// New scenario
				++scenNum;
				nFactors = iFactor;
				iFactor = 0;
			}
			
			if (buf.find("msmc 2000")!=string::npos ||
					buf.find("green")!=string::npos ||
					buf.find("Constant")!=string::npos ||
					buf.find("Trigger")!=string::npos ||
					buf.find("Term")!=string::npos) {  // Encontered new factor
				scenRecord1 = buf;
				
				using boost::is_any_of;                                                                                                                                     
				boost::algorithm::split(values, buf, is_any_of(",")); 
				curName = values[5];
				if(scenNum > 0)
					break;
				if(scenNum == 0) {
					curveName.push_back(curName);
					scenTerm1.push_back(values[13]);
				}
			}
			else {
				boost::algorithm::split(values, buf, is_any_of(",")); 
				// Only reading the scenarios for factors that either don't have term structure of have tenor 184 (for spot and vols)
				if (readEqOnly) {
					if(values.size() == 15 && (atoi(values[13].c_str())==0 || atoi(values[13].c_str())==184)) {
						
						scen_mat(scenNum, iFactor) = atof(values[14].c_str());
						++iFactor;
					}
				}
				else if(values.size() == 15){
						// change !!!
						if(scenNum == 0) {
							curNameTerm = curName + "_" + values[13];
							curveName.push_back(curNameTerm);
							scenTerm1.push_back(values[13]);
						}
						scen_mat(scenNum, iFactor) = atof(values[14].c_str());
						
						++iFactor;
					
				}
			}
		}
		else // end new scenario
			scenRecord1 = buf;
	} // end while


	// Statistics
	inp_file.close();
	scenNum += 1;
	out_file.open(out_fn);

	for(iFactor=0; iFactor < nFactors; ++iFactor) {
		out_file << curveName[iFactor] << ",";
	}
	out_file << endl;
	for(int i=0; i< scenNum; ++i) {
		for(iFactor=0; iFactor < nFactors; ++iFactor) {
				out_file << scen_mat(i,iFactor) << ",";
		}
		out_file << endl;
	}
	
	out_file.close();
	return 0;
}
*/
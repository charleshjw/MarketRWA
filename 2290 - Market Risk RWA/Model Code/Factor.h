#ifndef FACTOR_H
#define FACTOR_H

#define BOOST_DATE_TIME_NO_LIB // because of the Linker problem "cannot open file libboost_date_time-vc ... Before include headers
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
#include <numeric>      // std::inner_product

//#include <boost/date_time/gregorian/greg_month.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
//#include <boost/shared_ptr.hpp>
//#include <boost/range/numeric.hpp>  // for inner_product

class FactSet;
#include "ParamControl.h"
//#include "portfolio.h"
//#include "utilities.h"
#include <random>
using namespace boost::gregorian;
using namespace boost::numeric::ublas;

//std::vector<std::string> factTypeNames;

enum FactType {UnIdentified, FX_Forward, FX_Spot, FX_Pair_Vol, EQ_Vol, factor_GB, IRCap_Vol, IRSwaption_Vol, CDS, Turnover, Refi, MBS, IR_Curve, Eqt_Spot, num_FactType};
enum FactMappingType {NotMapped, Idio, Regressed, Proxied, PCA, NoHistData_NoProxy, NoHistData_ProxyOfDiffAssetClass, Unknown};
enum HistDataAvailabilityType {NoGaps, WithGaps, FewDates, NoData, Unidentified};

class CCAR_info;
class CCAR_shift;

class Factor {
public:
	
	Factor(int importance = 0) : _importance(importance), _spot(0.0), _resid_vol(0.0), _R2(0.0), _theTermForMapping(0), _var(0.0), _vol_from_mappedto_facs(0.0), _isLogNormal(true), _proxy_name("") {};
	Factor(std::string s_name); 
	Factor(int importance, unsigned int nTerms, FactType type, const std::vector<double>& vals, const std::vector<unsigned int>& terms, std::string curr, std::string subType, std::string s_name, std::string h_name) : 
	_importance(importance), _nTerms(nTerms), 	_type(type), _spot(vals), _terms(terms), _curr(curr), _subType(subType), _s_name(s_name), _h_name(h_name), _resid_vol(0), _R2(0), _theTermForMapping(0),  
		_mapping_type(Unknown), _vol_from_mappedto_facs(0.0), _isLogNormal(true), _proxy_name("") {};
	
	void			setNterms(unsigned int nTerms) {_nTerms = nTerms;};
	void			setType(FactType type) {_type = type;};
	void			setTerms(const std::vector<unsigned int>& terms) {_terms = terms; _nTerms = terms.size();};
	void			setCurr(std::string curr){_curr = curr;};
	void			setSubType(std::string subType) {_subType = subType; };
	void			setScenName(std::string name) {_s_name = name; };
	void			setHistName(std::string name) {_h_name = name; };
	//void setHistData(matrix<double>* matr) {_histData.reset(matr); };
	//void setRets(matrix<double>* rets) {_rets.reset(rets); };
	void			setHrDates(const std::vector<date>& hrDates) {_hrDates = hrDates;};
	void			setAttributes(std::map<std::string, std::string>& cur_2_reg_map, std::string hr_data_type, std::string factorGroup, std::string fullName, std::string curr, std::vector<date>& hr_dates, 
		std::map<date, std::vector<double> >& prc, std::map<date, std::vector<unsigned int> >& terms, int use_2_map_2 = 0);
	void			setAttributes(FactType iType, std::string subType, std::string cur, std::string reg, std::vector<double>& sens, std::vector<unsigned int> terms, int importance = 0);
	std::string		h_name_2_s_name_by_type(int type);
	std::string		h_name_2_s_name();
	void			setSpot(const std::vector<double>& prc);
	void			setForecast();
	void			setDelta(std::vector<double>& vals);
	void			setGamma(std::vector<double>& vals);
	void			setCrossGamma(matrix<double>& gam);
	void			setImportance(int importance) {_importance = importance;}; 
	void			setRegion(std::string reg) {_reg = reg;};
	bool			setRetsVols(ParamControl& params, std::vector<date>& dts);
	void			setVol(std::vector<double>& vol) {_vol = vol;};
	void			setR2(double R2) {_R2 = R2;};
	void			setRegrCoef(matrix<double>&  regr_coef) {_regr_coef = regr_coef;};
	void			setResidVol(double resid_vol) {_resid_vol = resid_vol;};
	date			setEODVal();
	void			assignEODVal(std::map<double, double>& mp, const FactSet& pos);
	void			setEODVal(std::vector<double>& vec) {_eod_val = vec;};
	void			setVar(double val) {_var = val;	};
	bool			setRetsOnDates(std::vector<date>& dts, unsigned int sim_horizon, matrix<double>& rets, std::vector<date>& new_dts);
	void			setMappingType(FactMappingType map_type) {_mapping_type = map_type;};
	void			setDelta(unsigned int iterm, double val);
	void			setProxyName(std::string h_name) { _proxy_name = h_name;};
	void			setProxyScalers(std::vector<double>& vec) {_proxy_scalers = vec;};
	void			setDataAvalType(HistDataAvailabilityType type) {_histDataAvalType = type;};
	void			setDistFlags(std::vector<bool>& flags) { _dist_flags = flags;};
	void			setUseProxyEOD(bool flag) {_useProxyEOD = flag;};
	void			setExcludeFromScenarios(bool val) {_excludeFromScenarios = val;};
	void			setIsLognormal(bool flag) {_isLogNormal = flag;};
	void			complete(std::vector<date>& dts);
	void			outputRegrCoef(std::ofstream& fn);
	unsigned int						getNterms() {
		return _nTerms;
	};
	std::vector<date>					getHrDates() {return _hrDates;};
	FactType							getType() {return _type;};
	std::string							getSubType() { return _subType;};
	const std::vector<unsigned int>&	getTerms() {return _terms;};
	std::string							getCurr() {return _curr;};
	std::string							getScenName() {return _s_name;};
	std::string							getHistName() {return _h_name;};
	std::string							getReg() {return _reg;};
	void								getTermsAtDate(date dd, std::vector<unsigned int>& terms);
	void								getDataAtDate(date dd, std::vector<double>& vals);
	std::vector<double>&				getVol() {return _vol;};
	std::vector<double>&				getDelta() {return _delta;};
	matrix<double>&						getGamma() {return _gamma;};
	matrix<double>&						getRegrCoef() {return _regr_coef;};
	double								getResidVol() {return _resid_vol;};
	double								getR2() {return _R2;};
	std::vector<boost::shared_ptr<Factor> >&	getMappedTo() {return _mapped_fact;};
	unsigned int						getTheTermForMapping() {return _theTermForMapping;};
	std::vector<double>&				getEODVal() {return _eod_val;};
	double								getVar() {return _var;};
	double								getMappedToVol() {return _vol_from_mappedto_facs;};
	bool								isLogNormal() {return _isLogNormal;};
	std::string							getProxyName() {return _proxy_name;};
	std::vector<unsigned int>&			getProxyTenors() {return _proxy_tenors;};
	std::vector<double>&				getProxyScalers() {return _proxy_scalers;};
	HistDataAvailabilityType			getDataAvalType() {return _histDataAvalType;};
	FactMappingType						getMappingType() {return _mapping_type;};
	void								genSims(matrix_range<matrix<double> >& rand_norm, std::vector<double>& wts);
	std::vector<bool>					getUseEmpirDist() {return _dist_flags;};
	bool								getUseProxyEOD() {return _useProxyEOD;};
	bool								getExcludeFromScenarios() {return _excludeFromScenarios;};

	bool								dataCleanUp(ParamControl& p);
	void								print_hist(std::ofstream& ff);
	void								print_factor_info(std::ofstream& ff);
	void								print_sim_rets(std::ofstream& ff);
	void								print_sim_diff(std::ofstream& ff);
	void								print_rets(std::ofstream& ff);
	void								print_vars(std::stringstream& ff);
	void								print_hist_pca(std::ofstream& ff);
	void								print_sensitivities(std::stringstream& ff);

	//boost::shared_ptr<matrix<double> >	getHistData();
	//boost::shared_ptr<matrix<double> >	getRets();
	void								factPrint(std::ofstream& of);
	int									comp_Scen_Hist_names(std::string curr, std::string subType, FactType type, std::string h_name);
	bool								calcPnL_HR(std::vector<double>& pnl);
	bool								calcPnL_MC(std::vector<double>& pnl, unsigned int var_percentile_int);
	void								calcAdhocShiftTilt(boost::numeric::ublas::vector<double>& shift, boost::numeric::ublas::vector<double>& tilt);
	bool								mapIt(ParamControl& params, std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Spot, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Spot, 
				std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_IR_Curve,
				std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Forward, 	std::vector<date>& dts, std::map<std::string, std::vector<std::string> >& reg_2_cur_map);
	bool								regress(std::vector<boost::shared_ptr<Factor> >& fvec, unsigned int sim_horizon, double decay);
	std::vector<double>&				getSpot() {return _spot;};
	//std::vector<boost::shared_ptr<Factor> >::iterator findProxyFact(std::vector<boost::shared_ptr<Factor> >& facs);
	double								interpolate(double x, std::string which_attribute);
	void								handle_noData_noProxy(double defaultIdioVol, double defaultIdioEquityVol);
	void								setRegBasedOnCur(const std::map<std::string, std::string>& cur_2_reg_map);
	// Related to CCAR
	void								scale_eod(std::map<double,double>& scalers);
	void								scale_rets(std::map<double, double>& scalers);
	void								save_eod();
	void								save_rets();
	void								reset_eod();
	void								reset_rets();
	void								apply_CCAR_info(const std::vector<boost::shared_ptr<CCAR_shift> >& ccar_shifts, unsigned int i_proj, const	std::map<std::pair<std::string, std::string>, 
		std::string>& mp, const std::map<std::string, std::string>& reg_2_cur_mp, const std::map<std::string, std::string>& cur_2_reg_map);
	//void								apply_CCAR_proj(const CCAR_info& info, unsigned int i_proj, std::map<std::string, std::string>& cur_2_reg_map);

	std::map<date, std::vector<double> >			_hrDataMap;
	std::map<date, std::vector<unsigned int> >		_hrTermsMap;
	matrix<double> 									_rets;
	matrix<double>									_rets_saved;
	matrix<double> 									_sim_rets;
	matrix<double> 									_histData;
	std::vector<unsigned int>						_mapped_terms;
	std::vector<double>								_delta;
	matrix<double>									_gamma;
	virtual ~Factor() {};
private:
	int										_importance;
	unsigned int							_nTerms;
	FactType								_type;
	std::vector<double>						_spot;
	std::vector<unsigned int>				_terms;
	std::string								_curr;
	std::string								_subType;
	std::string								_s_name;
	std::string								_h_name;
	std::string								_reg;
	std::vector<date>						_hrDates;
	
	std::vector<double>						_vol;
	std::vector<double>						_vol_saved;
	matrix<double>							_regr_coef; // has one column
	double									_resid_vol;
	double									_R2;
	std::vector<boost::shared_ptr<Factor> > _mapped_fact;
	unsigned int							_theTermForMapping;
	std::vector<double>						_eod_val;
	std::vector<double>						_eod_saved;
	double									_var;
	std::vector<double>						_var_by_term;
	FactMappingType							_mapping_type;
	double									_vol_from_mappedto_facs;
	bool									_isLogNormal; // if "false" is normal dynamics, otherwise assume log-normal
	std::string								_proxy_name;  // s_name
	std::vector<unsigned int>				_proxy_tenors;
	std::vector<double>						_proxy_scalers;
	HistDataAvailabilityType				_histDataAvalType;
	std::vector<bool>						_dist_flags;
	std::vector<double>						_vol_scalers_ccar;
	std::vector<double>						_eod_scalers_ccar;
	bool									_useProxyEOD;
	bool									_excludeFromScenarios;
	//bool									_isCoreRf;
	//bool									_isMultCrf;
	//matrix<double>						_ccar_shifts;
};


#endif FACTOR_H
#include "stdafx.h"

CoreFactor::CoreFactor(const Factor& ff) : Factor(ff) {
	//_hrDataMap = ff._hrDataMap;
}
/*
std::map<date, std::vector<double> >			_hrDataMap;
	std::map<date, std::vector<unsigned int> >		_hrTermsMap;
	matrix<double> 									_rets;
	matrix<double>									_rets_saved;
	matrix<double> 									_sim_rets;
	matrix<double> 									_histData;
	std::vector<unsigned int>						_mapped_terms;

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
	std::vector<double>						_delta;
	matrix<double>							_gamma;
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
	*/
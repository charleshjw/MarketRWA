#ifndef PORTFOLIO_H
#define PORTFOLIO_H

#include "Factor.h"
#include "CCAR_info.h"

enum ValFileType {Delta, CrossGamma};
enum ShiftType {Base, DeltaDn,  DeltaUp, CrossGammaUpUp, CrossGammaUpDn, CrossGammaDnUp, CrossGammaDnDn, numValues};

class ShiftScenInfo;

class FactSet {
public:
	FactSet() :  _systematicVar(0.0), _idioVar(0.0), _totalVar(0.0) {};
	FactSet(Factor& fac);
	FactSet(const std::vector<boost::shared_ptr<Factor> >& facs);
	FactSet(std::string name) : _name(name), _systematicVar(0.0), _idioVar(0.0), _totalVar(0.0) {};
	FactSet(const FactSet& pfl) : _name(pfl._name), _systematicVar(pfl._systematicVar), 
		_idioVar(pfl._idioVar), _totalVar(pfl._totalVar) {};
	double		getSystematicVar() {return _systematicVar;};
	double		getIdioVar() {return _idioVar;};
	double		getTotalVar() {return _totalVar;};
	std::string getName() {return _name;};
	const std::vector<boost::shared_ptr<Factor> >& getFactors() const {return _pos;};
	void			removeEqtFactors();
	void			readFromScenValFile(std::string fn, std::string info_fn, ValFileType type, std::vector<boost::shared_ptr<Factor> >& facs);
	void			setSensitivities(std::vector<boost::shared_ptr<Factor> >& facs); // set deltas equal to 1, zero gammas
	void			setSensitivities(std::string sens_fn, std::vector<boost::shared_ptr<Factor> >& facs);
	//std::string		setDeltasGammasFromShiftVals(std::vector<boost::shared_ptr<Factor> >& factors, std::map<std::string, std::vector<double> >& pfl_vals_deltas, double shift);
	//std::string		setCrossGammasFromShiftVals(std::vector<boost::shared_ptr<Factor> >& factors, std::map<std::string, std::vector<double> >& pfl_vals, double shift);
	void			readEOD(std::string dir);
	void			readEOD_equity_fx(std::string dir);
	void			getFactorsByTypeAndDataAval(HistDataAvailabilityType dataAvalType, FactType factType, std::vector<boost::shared_ptr<Factor> >& f_vec);
	void			getFactorsByMappingType(FactMappingType type, std::vector<boost::shared_ptr<Factor> >& o_facs);
	void			getFactorsByType(FactType factType, std::vector<boost::shared_ptr<Factor> >& f_vec);

	void			addPositions(std::vector<boost::shared_ptr<Factor> >& facs);
	void			setEOD_ProxiedUseProxyEOD();
	void			map_2_regress(ParamControl& params, std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Spot, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Spot, 
			std::vector<boost::shared_ptr<Factor> >& f_PCA_Eqt_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Vol, std::vector<boost::shared_ptr<Factor> >& f_PCA_IR_Curve,
			std::vector<boost::shared_ptr<Factor> >& f_PCA_FX_Forward, 	std::vector<date>& dts, std::map<std::string, std::vector<std::string> >& reg_2_cur_map, std::vector<date>& dts_common);
	
	std::vector<boost::shared_ptr<Factor> >::const_iterator		findFactor(std::string name, bool isHistName) const;
	//std::vector<boost::shared_ptr<Factor> >::iterator		matchNameWithFactor(std::string name, bool isHistName);
	std::vector<boost::shared_ptr<Factor> >::const_iterator		findProxyFact(std::string proxy_name) const;

	void		genSimsProxied(std::vector<boost::shared_ptr<Factor> >& proxied, unsigned int nsim);

	void		assignEmpirDistFlag(std::map<std::string, std::vector<unsigned int> >& empirDistMap);
	void		setName(std::string name) {_name = name;};
	double		calc_MC_Var(double percentile);
	
	void		print_vars(std::string fn);
	void		print_sensitivities(std::string fn);
	
	void		output_factors_info(std::string dir_output);
	void		output_mapping_info(std::string dir_output);
	void		output_sim_diff(std::string fn);
	void		output_hist_rets(std::string fn);
	int			output_scenarios2(std::string dir, std::string asOfDate, unsigned int nsim, bool output_diff);
	void		output_hr_factors(std::string fn, unsigned int sim_horizon, double decay);
	void		calc_factor_coef(std::map<unsigned int, std::string>& stockMap, std::map<unsigned int, std::vector<double> >& factorCoefMap, 
					std::vector<boost::shared_ptr<Factor> >& pca, unsigned int sim_horizon, double decay);
	void		output_factor_coef(std::map<unsigned int, std::string>& stockMap, std::map<unsigned int, std::vector<double> >& factorCoefMap, std::string fn);

	void		readProxyMap(std::string fn, double defaultIdioVol, double defaultIdioEquityVol);
	void		save_rets_eod();
	void		reset_rets_eod();
	void		apply_CCAR_info(const std::vector<boost::shared_ptr<CCAR_shift> >& ccar_shifts, unsigned int i_proj, const std::map<std::pair<std::string, std::string>, std::string>& rfTypeCur_coreRf_mp, 
					const std::map<std::string, std::string>& reg_2_cur_mp, const std::map<std::string, std::string>& cur_2_reg_map);
	void		clear_CCAR_info(const CCAR_info& info, unsigned int ii);
	void		digest_crf_info(const CCAR_info& ccar_info);
	void		genSensShifts_all(std::string asOfDate, std::string dir_shifts, double shift, unsigned int npca);
	//int			genShifts_delta(std::string file_nm, std::string factors_info_fn, std::string asOfDate, double shft1bp, std::string scenName);
	int			genShifts(std::string file_nm, std::string factors_info_fn, std::string asOfDate, double shft1bp, std::string scenName, bool doIR, bool genShifts, bool getShiftsInfo);
	//int			genShifts_cross_gamma(FactType factType, std::string file_nm, std::string factors_info_fn, std::string asOfDate, double shift);
	void		output_scenario(std::string scenSetName, std::string scen_name, std::string scen_color, std::string ScenarioStartTime, std::string TimeEvolution, std::string TimeEvolutionToTrigger,
			std::string TriggerHolder, std::string ScenShiftRule, std::string  ScenType,std::string  ScenReplVal, std::string asOfDate, 
			ShiftType shiftType, unsigned int k, double prob, double shift, bool doIR, std::vector<boost::shared_ptr<Factor> >::iterator it, 
			std::vector<boost::shared_ptr<Factor> >::iterator jt, unsigned int iterm, unsigned int jterm, std::stringstream & out_file);
	bool		setSens_fromShiftVals(const std::map<unsigned int, ShiftScenInfo>& mp, std::map<unsigned int, double>& val_mp, double shift);
	void		resetSens();
	void		makeSinglePCAFactor(unsigned int num_pca);
	void		readSims(std::string fn, unsigned int nsim);

	//void		set_and_eliminate_factType(const FactSet & pfl, FactType factType);
	~FactSet() {};

private:
	std::string _name;
	
	double _systematicVar;
	double _idioVar;
	double _totalVar;
	std::vector<boost::shared_ptr<Factor> > _pos;
};

#endif PORTFOLIO_H
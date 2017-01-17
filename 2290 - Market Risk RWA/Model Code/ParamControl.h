#ifndef PARAMCONTROL_H
#define PARAMCONTROL_H

#define BOOST_DATE_TIME_NO_LIB // because of the Linker problem "cannot open file libboost_date_time-vc ... Before include headers

#include <string>
#include <iostream>
#include <fstream>
#include<boost/algorithm/string/split.hpp>                                      
#include<boost/algorithm/string.hpp>

class ParamControl {
public:
	ParamControl () : _minDatesCount(450), _num_pca_ir_curve(18), _num_pca_fx_forward(18), _num_pca_eqt_spot(18), _num_pca_eqt_vol(18), _num_pca_fx_spot(18), _num_pca_fx_vol(18), _num_pca_ircap_vol(3),
		_minNumRetsForRegr(21), _maxAllowedGap(5), _maxDatesUnch(10), _defaultIdioVol(0.05), _defaultIdioEquityVol(0.05), _minPositiveValue(1E-4), _staleDataDailyVolMin(0.0),
	_unreasonableValue(-123456), _maxDates(3600), _maxTerms(100), _factorCoefMax(10000), _MC_num(1000), _decay_factor(1.0), _percentile(0.01), _randomness_control(0),
	_largestAllowedReturn(0.5), _shift(0.0001), _sim_horizon(1) {};
	~ParamControl () {};
	void set_minDatesCount(unsigned int v) {_minDatesCount = v;};
	void set_num_pca(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4, unsigned int n5, unsigned int n6) {
		_num_pca_ir_curve = n1; 
		_num_pca_fx_forward = n2;
		_num_pca_eqt_spot = n3;
		_num_pca_eqt_vol = n4;
		_num_pca_fx_spot = n5;
		_num_pca_fx_vol = n6;};
	void set_minNumDatesForRegr(unsigned int v) {_minNumRetsForRegr = v;};
	void set_maxAllowedGap(unsigned int v) {_maxAllowedGap = v;};
	void set_maxDatesUnch(unsigned int v) {_maxDatesUnch = v;};
	void set_defaultIdioVol(double v) {_defaultIdioVol = v;};
	void set_defaultIdioEquityVol(double v) {_defaultIdioEquityVol = v;};
	void set_minPositiveValue(double v) {_minPositiveValue = v;};
	void set_staleDataDailyVolMin(double v) {_staleDataDailyVolMin = v;};
	void set_unreasonableValue(double v) {_unreasonableValue = v;};
	void set_maxDates(unsigned int v) {_maxDates = v;};
	void set_maxTerms(unsigned int v) {_maxTerms = v;};
	void set_factorCoefMax(unsigned int v) {_factorCoefMax = v;};
	void setNsim(unsigned int nsim) {_MC_num = nsim;};

	unsigned int get_minDatesCount() {return _minDatesCount;};
	unsigned int get_minNumRetsForRegr() {return _minNumRetsForRegr;};
	unsigned int get_maxAllowedGap() {return _maxAllowedGap;};
	unsigned int get_maxDatesUnch() {return _maxDatesUnch;};
	double		get_defaultIdioVol() {return _defaultIdioVol; };
	double		get_defaultIdioEquityVol() {return _defaultIdioEquityVol; };
	double		get_minPositiveValue() {return _minPositiveValue;};
	double		get_staleDataDailyVolMin() {return _staleDataDailyVolMin;};
	double		get_unreasonableValue() {return _unreasonableValue;};
	unsigned int get_maxDates() {return _maxDates;};
	unsigned int set_maxTerms() {return _maxTerms;};
	unsigned int get_factorCoefMax() {return _factorCoefMax;};
	unsigned int getNsim() {return _MC_num;};
	unsigned int get_npca_ir_curve() {return _num_pca_ir_curve;};
	unsigned int get_npca_fx_forward() {return _num_pca_fx_forward;};
	unsigned int get_npca_eqt_spot() {return _num_pca_eqt_spot;};
	unsigned int get_npca_eqt_vol() { return _num_pca_eqt_vol;};
	unsigned int get_npca_fx_spot() {return _num_pca_fx_spot;};
	unsigned int get_npca_fx_vol() { return _num_pca_fx_vol;};
	unsigned int get_npca_itcap_vol() { return _num_pca_ircap_vol;};
	double		 get_decay_factor() { return _decay_factor;};
	double		 get_percentile() { return _percentile;};
	unsigned int get_randomness_control() { return _randomness_control;};
	double		get_maxAllowedRet() { return _largestAllowedReturn;};
	double		getShiftForGreeks() { return _shift;};
	unsigned int getSimHorizon() {return _sim_horizon;};
	void readFromFile(std::string fnm);

private:
	unsigned int	_minDatesCount;			// Minimum number of dates a factor should have to avoid being regressed
	unsigned int	_num_pca_ir_curve;		// Number of principla components calculated for equities to use in regressions
	unsigned int	_num_pca_fx_forward;
	unsigned int	_num_pca_eqt_spot;
	unsigned int	_num_pca_eqt_vol;
	unsigned int	_num_pca_fx_spot;
	unsigned int	_num_pca_fx_vol;
	unsigned int	_num_pca_ircap_vol;
	unsigned int	_minNumRetsForRegr;	// Minimum days of data without gap a factor should have to be regressed
	unsigned int	_maxAllowedGap;			// Maximum allowed gap in the data. If a gap exceeds this number, all previous values are thrown away
	unsigned int	_maxDatesUnch;			// Maximum allowed number of dates with unchanged data. If a number of days with unchanged data exceeds this number, all previous values are thrown away
	double			_defaultIdioVol;
	double			_defaultIdioEquityVol;
	double			_minPositiveValue;
	double			_staleDataDailyVolMin;
	double			_unreasonableValue;
	unsigned int	_maxDates;
	unsigned int	_maxTerms;
	unsigned int	_factorCoefMax;
	unsigned int	_MC_num;				// number of MC simulations
	double			_decay_factor;
	double			_percentile;
	unsigned int	_randomness_control;
	double			_largestAllowedReturn;	// Returns larger than this are considered errouneous
	double			_shift;					// Shifts to set scenarios for sensitivities
	unsigned int	_sim_horizon;			// simulation horizon in days
};

#endif PARAMCONTROL_H
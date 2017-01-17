// NewScenGen.cpp : Defines the entry point for the console application.
//

// aj - include the header file that has all the header files
#include "stdafx.h"

#ifdef _DEBUG
#define DEBUG_NEW new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define new DEBUG_NEW
#endif

using namespace boost;
using namespace boost::assign; // bring 'operator+=()' into scope
using namespace boost::numeric::ublas;
using namespace std;

extern std::vector<std::string> factTypeNames;

string	readShiftScen_Values(std::string fn, std::map<std::string, std::vector<double> >& pfl_val_map);
string	readShiftScen_FactorInfo(std::string info_fn, std::vector<std::string>& fact_names, std::vector<unsigned int>& nterms, double& shift);

void	imposeCorr_REFI_TURNOVER(std::vector<boost::shared_ptr<Factor> >& f_MBS_Adhoc, double corr);
string	set_MBS_Ad_Hoc_RF_params(string dir_hist_data, std::vector<boost::shared_ptr<Factor> >& facs, double& cor);

int		read_scen(std::string  f_name, int nscen, std::vector<string>& names, matrix<double>&mat);
int		readScenarioFile(string inp_fn, std::vector<boost::shared_ptr<Factor> >& factors, unsigned int nsim);

double	getVarFromPnl(std::vector<double>& pnl, double perc);
void	populateFactor(Factor* factor, std::vector<string> nameParts, int& ic1, int& ic2, int& ic3, int& ic4, int& ic6, int& unIdentified);
int		readHistData(ErrorLog& errLog, ParamControl& params, std::map<string, string>& cur_2_reg_map, string dir, std::vector<string>& h_names, std::vector<boost::shared_ptr<Factor> >& factors);
int		readHistData_clean(ErrorLog& errLog, ParamControl& params, std::map<string, string>& cur_2_reg_map, string dir, std::vector<string>& h_names, std::vector<boost::shared_ptr<Factor> >& factors);
void	assignCCAR_Mapping(std::vector<boost::shared_ptr<Factor> >& factors);

int		alignDates(ErrorLog& errLog, ParamControl& params, std::vector<boost::shared_ptr<Factor> >& factors, std::vector<date>& dts_common);
int		alignSens_HistFactors(std::vector<boost::shared_ptr<Factor> >& s_factors, std::vector<boost::shared_ptr<Factor> >& h_factors);
double	calc_Hist_Var(std::vector<boost::shared_ptr<Factor> >& factors, std::vector<date>& dts, double percentile);
int		genShifts2CalcSens(string inp_scen_fn, string shifts_dir, string factor_info_fn);

double doAll(string dir_do_all, string asOfDate, bool do_GenSims, bool do_GenShifts, bool do_CalcVar, bool readFactSetValuations);

unsigned int getNumFactors(std::vector<boost::shared_ptr<Factor> >& f_vec);

// Maps
int read_cur_2_reg_map(string map_fn, std::map<string, string>& cur_2_reg_map);
int invert_cur_2_reg(std::map<string, string>& cur_2_reg_map, std::map<string, std::vector<string> >& reg_2_cur_map);

double genSimsOld(std::vector<boost::shared_ptr<Factor> >& f_notMapped, std::vector<boost::shared_ptr<Factor> >& f_Regressed, std::vector<boost::shared_ptr<Factor> >& f_Idio, 
	std::vector<boost::shared_ptr<Factor> >& f_PCA, std::vector<date>& dts_common, unsigned int nsim, std::vector<double>& wts, unsigned int sim_horizon = 1, unsigned int randomness_control = 0);
void	genSimsProxied(std::vector<boost::shared_ptr<Factor> >& proxied, std::vector<boost::shared_ptr<Factor> >& facs, unsigned int nsim);

double calc_pca_projections(std::vector<boost::shared_ptr<Factor> >& facs, unsigned int n_eig, double decay, std::vector<boost::shared_ptr<Factor> >& f_PCA, 
	std::vector<boost::shared_ptr<Factor> >& f_PCA_All, unsigned int sim_horizon);
void calc_hr_eqt_factors_and_coef(std::vector<boost::shared_ptr<Factor> >& f_Eq_Spot, double decay, unsigned int sim_horizon, std::string fn_hr, 
	std::string fn_coef, std::vector<boost::shared_ptr<Factor> >& f_Eq_Spot_Mapped);
void read_pfl_list(std::string fn, std::vector<std::string>& list);
void output_valuations(std::string fn, const std::map<std::string, std::map<unsigned int, double> >& mp_all);

void read_RF_names(string dir, std::vector<string>& h_names);
void readShiftInfo(std::string str, std::map<unsigned int, ShiftScenInfo>& mp);
bool readVals(std::string str, std::map<unsigned int, double>& val_mp);
bool readVals_all(std::string str, std::map<std::string, std::map<unsigned int, double> >& val_mp_all, std::vector<std::string>& list);
string eqt_idx_h_names_lst[] = {"C12496G10_OTC_RUSSELL3000Spot"};
// "C85299600_OTC_SPTRSpot", "S1782608_NYSE_GDDUEAFESpot", "S1788974_NYSE_GDDUWISpot","S1800942_NYSE_NDUEACWFSpot", "S1782620_NYSE_NDDUEAFESpot", "S1800793_NYSE_NDUEEGFSpot"
static const std::vector<string> eqt_idx_h_names(eqt_idx_h_names_lst, eqt_idx_h_names_lst + sizeof(eqt_idx_h_names_lst)/sizeof(eqt_idx_h_names_lst[0]));

std::ofstream eig_info, log_ff; // temp Tanya

//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, char* argv[])
{
	string dir_do_all, asOfDate;
	bool do_GenSims(false), do_GenShifts(false), do_CalcVar(false), readFactSetValuations(false), do_CCAR(false), do_removeEqtFactFromScenarios(false), do_output_MC_scen(false), do_ReadSims(false);

	dir_do_all = std::string(argv[1]);
	asOfDate = std::string(argv[2]);

	if((*argv[3]) == 'Y')
		do_GenSims = true;
	if((*argv[4]) == 'Y')
		do_removeEqtFactFromScenarios = true;
	if((*argv[5]) == 'Y')
		do_GenShifts = true;
	if((*argv[6]) == 'Y')
		do_CalcVar = true;
	if((*argv[7]) == 'Y')
		readFactSetValuations = true;
	if((*argv[8]) == 'Y')
		do_CCAR = true;
	if((*argv[9]) == 'Y')
		do_output_MC_scen = true;
	if((*argv[10]) == 'Y')
		do_ReadSims = true;

	FactSet pfl;

	std::vector<std::map<unsigned int, ShiftScenInfo> > shift_info_mp;
	std::vector<std::map<std::string, std::map<unsigned int, double> > > val_mp_all; // prortfolio name is mapped to valuation results for all sensitivity shift scenarios
	std::map<std::string, std::map<unsigned int, double> >::iterator itmap1, itmap2;
	std::string pfl_val_names[] = {"val_non_ir.csv", "val_ir.csv"}, shift_info_names[] = {"shifts_non_ir_info.csv", "shifts_ir_info.csv"};
	std::vector<std::string>	pfl_val(std::vector<std::string>(pfl_val_names, pfl_val_names + sizeof(pfl_val_names)/sizeof(std::string))), 
							 shift_info(std::vector<std::string>(shift_info_names, shift_info_names + sizeof(shift_info_names)/sizeof(std::string)));
	val_mp_all.resize(pfl_val.size());
	shift_info_mp.resize(pfl_val.size());

	double var_MC(1234), var_HR(4321), defaultIdioVol(0.05), defaultIdioEquityVol(0.05), shift(0.0001), percentile(0.01);
	int r1(0), r2(0), r4(0), r5(0), m1(0), m2(0), nrets(0);
	unsigned int minNumRetsForRegress, nSim, i, randomness_control(0), sim_horizon(1), ii, n_proj;
	bool histVar(false), cleanExists(false), res(false);;
	date eod_date(1900,1,1);
	std::map<string, string> cur_2_reg_map;
	std::map<string, std::vector<string> > reg_2_cur_map;
	std::vector<boost::shared_ptr<Factor> > s_factors;	// factors from scenario file
	std::vector<boost::shared_ptr<Factor> > h_factors, f_Eqt_Index, f_Eq_Spot, f_Eq_Spot_Mapped, f_EQ_Vol, f_FX_Vol, f_IRCap_Vol, 
		f_PCA_Eqt_Spot, f_PCA_FX_Spot, f_PCA_Eqt_Vol, f_PCA_FX_Vol, f_PCA_IR_Curve, f_PCA_FX_Forward, f_PCA_IRCap_Vol, f_MBS_Adhoc;	// factors from hist data
	std::vector<boost::shared_ptr<Factor> > f_notMapped, f_Regressed, f_Idio, f_Proxied, f_Not_Proxied;;
	std::vector<boost::shared_ptr<Factor> > f_PCA;
	std::vector<std::string> company_tripples_list;

	std::vector<date> dts_common;
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::vector<std::string> h_names;
	std::ofstream scen_file, cleaned_facs, param_controle_file, log_file;
	
	std::vector<double> wts;
	std::vector<std::string> empir_dist_fact_names;
	std::vector<boost::shared_ptr<CCAR_shift> > ccar_shifts;
	CCAR_info ccar_info;
	FactSet pfl_PCA, pfl_CCAR;

	// Set Input names
	string dir_rf(dir_do_all + "RF\\");					// List yesterday's positions risk factors
	string dir_hist_data(dir_do_all + "HistData\\");	// Historical data
	string dir_input(dir_do_all + "Input\\");			// Param control file
					
	// Set output directories names
	string dir_output(dir_do_all + "Output\\");
	string dir_shifts(dir_do_all + "Shifts\\");
	string dir_maps(dir_do_all + "Maps\\");
	string dir_logs(dir_do_all +"Logs\\");
	string dir_sens(dir_do_all +"Sens\\");
	string dir_var(dir_do_all +"Var\\");
	string dir_scenarios(dir_do_all +"Scenarios\\");
	string dir_pca(dir_do_all +"PCA\\");
	string dir_pfl_vals(dir_do_all +"Pfl_vals\\");
	string dir_proxy_maps(dir_do_all + "ProxyMaps\\");
	string dir_eod_data(dir_do_all + "EODdata\\");
	string dir_CCAR(dir_do_all + "CCAR\\");
	string cleaned_facs_fn(dir_output + "cleaned_facs.csv"), control_params_fn(dir_input +"control_params.csv"), maps_fn(dir_maps + "cur_2_reg_map.csv"), stockmap_fn(dir_input +"stockmap.dat");
	string pfl_name("");

	std::string sens_file_name("sens_test.csv");
	log_ff.open(dir_output + "log_test.csv");

	std::ofstream test_f;
	test_f.open(dir_output + "test.csv");

	// read stock map file
	std::map<unsigned int, std::string> stockMap;
	std::map<unsigned int, std::vector<double> > factorCoefMap;
	readStockmapFile(stockmap_fn, stockMap);

	// Read control parameters
	ParamControl params;
	params.readFromFile(control_params_fn);
	defaultIdioVol =		params.get_defaultIdioVol();
	defaultIdioEquityVol =	params.get_defaultIdioEquityVol();
	minNumRetsForRegress =	params.get_minNumRetsForRegr();
	nSim =					params.getNsim();
	shift =					params.getShiftForGreeks();
	sim_horizon =			params.getSimHorizon();

	// If read scenario file and convert to matrix form
	if(1) {
		
		/*
		std::vector<boost::shared_ptr<Factor> > facs;
		readScenarioFile(dir_scenarios + "varscen_ori.csv", facs, nSim + 1);
		//readScenarioFile(dir_scenarios + "scen_MC.csv", facs, nSim + 1);
		ofstream ff;
		//ff.open(dir_scenarios + "varscen_ori_matrix.csv");
		ff.open(dir_scenarios + "scen_MC_matrix.csv");
		
		for(itfac = facs.begin(); itfac != facs.end(); ++itfac) 
			(*itfac)->print_sim_rets(ff);
		ff.close();
		
		read_pfl_list(dir_do_all + "COMPANY_reg_tripples_list.csv", company_tripples_list);
		res = readVals_all(dir_pfl_vals + pfl_val[1], val_mp_all[1], company_tripples_list);	
		output_valuations(dir_pfl_vals + "clean_" + pfl_val[1], val_mp_all[1]);
		*/
		
		read_pfl_list(dir_do_all + "COMPANY_reg_tripples_list.csv", company_tripples_list);
			
			for(ii = 0; ii < pfl_val.size(); ++ii) {

				// Read valuation results and shifts info
				cleanExists = readVals_all(dir_pfl_vals + "clean_" + pfl_val[ii], val_mp_all[ii], company_tripples_list);
				if(!cleanExists) {	
					res = readVals_all(dir_pfl_vals + pfl_val[ii], val_mp_all[ii], company_tripples_list);
					if(!res)	
						return 9991;
					output_valuations(dir_pfl_vals + "clean_" + pfl_val[ii], val_mp_all[ii]);
				}
				//readShiftInfo(dir_shifts + shift_info[ii], shift_info_mp[ii]);
			}
		
		return 9999;
	}

	// Read currency to region map - these are not used actually 
	m1 = read_cur_2_reg_map(maps_fn, cur_2_reg_map);		// each currency is roughly mapped to its geographical region
	m2 = invert_cur_2_reg(cur_2_reg_map, reg_2_cur_map);

	if(!do_ReadSims) {
		// Read file with the list of factors for which Empirical Dist will be used
		std::string fn_empir_dist(dir_input + "realFactor.csv");
		std::map<std::string, std::vector<unsigned int> > empirDistMap;
		readEmpirDistFactList(fn_empir_dist, empirDistMap);

		// Read previous day risk factor names to identify relevant factors
		read_RF_names(dir_rf, h_names);

		// Read historical data, construct factors
		ErrorLog errLog_readHistData(dir_logs + "err_readHistData.csv");
		//r2 = readHistData(errLog_readHistData, params, cur_2_reg_map, dir_hist_data, h_names, h_factors);
		r2 = readHistData_clean(errLog_readHistData, params, cur_2_reg_map, dir_hist_data, h_names, h_factors);
		errLog_readHistData.print();

		// Align dates, create clean set to serve as risk factors used to map other risk factors to deal with factors that don't have sufficient data
		// Set HistDataAvalType to NoGaps and WithGaps
		// This routine might reduce the number of terms if the data is poor
		ErrorLog errLog_alignDates(dir_logs + "err_alignDates.csv");
		double decay(params.get_decay_factor());
		r4 = alignDates(errLog_alignDates, params, h_factors, dts_common);
		nrets = dts_common.size() - 1;

		for(itfac = h_factors.begin(); itfac != h_factors.end(); ++itfac)
			if((*itfac)->getNterms() == 0)
				int tt = 0.0;

		// Check if some risk factors from the list don't have historical data and include them as well
		if(1) { // tempTanya
			for(i = 0; i < h_names.size(); ++i) {
				itfac = h_factors.begin();
				while(itfac != h_factors.end() && h_names[i] != (*itfac)->getHistName())
					++itfac;
				if(itfac == h_factors.end()) { // Create a factor which will be treated as idiosyncratic
					if(h_names[i] == "C00287Y10_NYSE_ABBVSpot")
						int tt = 0; // temp Tanya

					boost::shared_ptr<Factor> ptr(new Factor);
					std::string s_name("");

					(*ptr).setNterms(0);
					(*ptr).setType(UnIdentified);
					(*ptr).setHistName(h_names[i]);
					s_name = (*ptr).h_name_2_s_name();
					(*ptr).setScenName(s_name);

					(*ptr).setCurr("USD");
					(*ptr).setRegion("North America");
					(*ptr).setR2(0.0);
					(*ptr).setResidVol(0.0);
					(*ptr).setProxyName("");
					(*ptr).setUseProxyEOD(false);
					(*ptr).setMappingType(Idio);
					(*ptr).setDataAvalType(NoData);
					(*ptr).setExcludeFromScenarios(false);

					if(h_names[i] == "REFI_DAILY" || h_names[i] == "TURNOVER_DAILY") {
						(*ptr).setType(MBS);
						f_MBS_Adhoc.push_back(ptr);
					}
					h_factors.push_back(std::move(ptr));
				}
			}
		}
		pfl.addPositions(h_factors);

		// Read Proxy Table
		std::map<std::string, std::pair<std::string, double> > proxy_map;
		std::string proxy_fn(dir_proxy_maps + "RWLink_Proxies.csv"); // For VaR
		pfl.readProxyMap(proxy_fn, defaultIdioVol, defaultIdioEquityVol);

		// Read EOD data
		pfl.readEOD(dir_eod_data);
		pfl.readEOD_equity_fx(dir_eod_data);

		// Handle Proxied that need to use Proxy EOD
		//pfl.setEOD_ProxiedUseProxyEOD();

		//if(do_GenSims && dts_common.size() >= params.get_minDatesCount()) {
		if(do_GenSims) {
			double pca_eq_pct(0), pca_fx_spot_pct(0.0), pca_fx_vol_pct(0.0), pca_ir(0.0), pca_fx_fwd(0.0), pca_ircap_vol(0.0);
			// Sort "good" factors for the following factor types FX_Forward, FX_Spot, IR_Curve, Eqt_Spot, Eqt_Vol and FX_Vol to be used for regression
			std::vector<boost::shared_ptr<Factor> > f_Eqt_Spot, f_IR_Curve, f_FX_Spot, f_FX_Forward;
			pfl.getFactorsByTypeAndDataAval(NoGaps, Eqt_Spot, f_Eq_Spot);
			pfl.getFactorsByTypeAndDataAval(NoGaps, IR_Curve, f_IR_Curve);
			pfl.getFactorsByTypeAndDataAval(NoGaps, FX_Spot, f_FX_Spot);
			pfl.getFactorsByTypeAndDataAval(NoGaps, FX_Forward, f_FX_Forward);
			//getFactorsByType(h_factors, EQ_Vol, f_EQ_Vol);
			pfl.getFactorsByTypeAndDataAval(NoGaps, FX_Pair_Vol, f_FX_Vol);
			pfl.getFactorsByTypeAndDataAval(NoGaps, IRCap_Vol, f_IRCap_Vol);

			// Create pca projections
			eig_info.open(dir_output + "eig_info.csv");
			pca_eq_pct = calc_pca_projections(f_Eq_Spot, params.get_npca_eqt_spot(), decay, f_PCA_Eqt_Spot, f_PCA, sim_horizon);
			pca_fx_spot_pct = calc_pca_projections(f_FX_Spot, params.get_npca_fx_spot(), decay, f_PCA_FX_Spot, f_PCA, sim_horizon);
			pca_fx_vol_pct = calc_pca_projections(f_FX_Vol, params.get_npca_fx_vol(), decay, f_PCA_FX_Vol, f_PCA, sim_horizon);
			pca_ir = calc_pca_projections(f_IR_Curve, params.get_npca_ir_curve(), decay, f_PCA_IR_Curve, f_PCA, sim_horizon);
			pca_fx_fwd = calc_pca_projections(f_FX_Forward, params.get_npca_fx_forward(), decay, f_PCA_FX_Forward, f_PCA, sim_horizon);
			pca_ircap_vol = calc_pca_projections(f_IRCap_Vol, f_IRCap_Vol.size(), decay, f_PCA_IRCap_Vol, f_PCA, sim_horizon);

			pfl.map_2_regress(params, f_PCA_Eqt_Spot, f_PCA_FX_Spot, f_PCA_Eqt_Vol, f_PCA_FX_Vol, f_PCA_IR_Curve, f_PCA_FX_Forward, 
				dts_common, reg_2_cur_map, dts_common);
		

			if(f_Eq_Spot.size() > 0) {
				pfl_PCA.setSensitivities(f_PCA_Eqt_Spot);
			
				std::string fn(dir_pca +"hr_factorGB.csv");
				pfl_PCA.output_hr_factors(fn, sim_horizon, decay);
			
				pfl.calc_factor_coef(stockMap, factorCoefMap, f_PCA_Eqt_Spot, sim_horizon, decay);
				eig_info << "f_Eq_Spot.size," <<f_Eq_Spot.size() << "," << stockMap.size() << "," << factorCoefMap.size() << std::endl;
				std::string factorCoef_fn(dir_pca + "factorCoefGB.csv");
				pfl.output_factor_coef(stockMap, factorCoefMap, factorCoef_fn);

			}
			eig_info.close();

			// Output historical factors info
			pfl.output_factors_info(dir_output);
		
			// Output mapping information for the mapped factors
			pfl.output_mapping_info(dir_output);

			// Generate scenarios
			calculate_weights(decay, dts_common, wts, sim_horizon);
			randomness_control = params.get_randomness_control();

			pfl.getFactorsByMappingType(NotMapped, f_notMapped);
			pfl.getFactorsByMappingType(Regressed, f_Regressed);
			pfl.getFactorsByMappingType(Idio, f_Idio);
			pfl.getFactorsByMappingType(Proxied, f_Proxied);
	
			// Add equity PCA to the set
			pfl.addPositions(f_PCA_Eqt_Spot);

			// Assign empirical distribution flag
			pfl.assignEmpirDistFlag(empirDistMap);
		
			genSimsOld(f_notMapped, f_Regressed, f_Idio, f_PCA, dts_common, nSim, wts, sim_horizon, randomness_control);

			// Handle REFI and TURNOVER - uncorrelated with everything else but correlated with each other
			double corr(-0.5);
			set_MBS_Ad_Hoc_RF_params(dir_hist_data, f_MBS_Adhoc, corr);
			imposeCorr_REFI_TURNOVER(f_MBS_Adhoc, corr);

			// Handle Proxied factors
			pfl.genSimsProxied(f_Proxied, nSim);
		
			// Output historical returns
			pfl.output_hist_rets(dir_output + "hist_rets.csv");

			// Outpot PCA
			//pfl_PCA.output_sim_diff(dir_pca + "PCA_EQ_Spot_sim.csv");

			if(do_removeEqtFactFromScenarios)
				pfl.removeEqtFactors();
		
			// Output simulation results
			pfl.output_sim_diff(dir_output + "scen_diff.csv");

			if(do_output_MC_scen)  
				pfl.output_scenarios2(dir_scenarios, asOfDate, nSim, true);
		
			if(0) 
				test_hist_vs_simulated(dir_output, f_notMapped);
			if(0)
				test_hist_vs_simulated_Mapped(dir_output, f_notMapped, f_Regressed, sim_horizon);
		}
		
	}
	else { // read sims
		pfl.readSims(dir_output + "scen_diff.csv", nSim);
		pfl.output_sim_diff(dir_output + "scen_diff_new.csv");
	}
	if(do_GenShifts || readFactSetValuations) {	
		// Get rid of "_Return" curves for equity
		pfl.removeEqtFactors();

		if(do_GenShifts) {
			pfl.genSensShifts_all(asOfDate, dir_shifts, shift, params.get_npca_eqt_spot());
			pfl.output_sim_diff(dir_output + "scen_diff2.csv");
		}
		if(readFactSetValuations) { // Read factor info file in Shifts directory and portfolio valuation results in Pfl_val directory to calculate and set sensitivities

			// Put PCA in one factor with multiple tenors
			pfl.makeSinglePCAFactor(params.get_npca_eqt_spot());
			pfl.output_sim_diff(dir_output + "scen_diff2.csv");

			// Read production triples
			read_pfl_list(dir_do_all + "COMPANY_reg_tripples_list.csv", company_tripples_list);
			
			for(ii = 0; ii < pfl_val.size(); ++ii) {

				// Read valuation results and shifts info
				cleanExists = readVals_all(dir_pfl_vals + "clean_" + pfl_val[ii], val_mp_all[ii], company_tripples_list);
				if(!cleanExists) {	
					res = readVals_all(dir_pfl_vals + pfl_val[ii], val_mp_all[ii], company_tripples_list);
					if(!res)	
						return 9991;
					output_valuations(dir_pfl_vals + "clean_" + pfl_val[ii], val_mp_all[ii]);
				}
				readShiftInfo(dir_shifts + shift_info[ii], shift_info_mp[ii]);
			}
			std::ofstream ff_var;
			ff_var.open(dir_var + "Company_VaR_ir.csv");
			// For triples
			i = 0;
			//itmap2 = val_mp_all[1].begin();

			for(itmap2 = val_mp_all[1].begin(); itmap1 != val_mp_all[1].end(); ++itmap2) {
					
				pfl.resetSens();
				pfl_name = itmap1->first;
				// Non-ir
				//pfl.setSens_fromShiftVals(shift_info_mp[0], itmap1->second, shift);
				// Ir
				pfl.setSens_fromShiftVals(shift_info_mp[1], itmap2->second, shift);
				
				std::ostringstream ostr;
				ostr << i;

				pfl.print_sensitivities(dir_output + "sensitivities" + ostr.str() + ".csv");
				var_MC = pfl.calc_MC_Var(percentile);
					
				ff_var << pfl_name << "," << var_MC << std::endl;
				pfl.print_vars(dir_var + "var_report_" + ostr.str() + ".csv");
				++i;
				//++itmap2;
			}
			// For the top level portfolio
			//itmap2 = val_mp_all[1].begin();
			pfl.resetSens();
			for(itmap2 = val_mp_all[1].begin(); itmap1 != val_mp_all[1].end(); ++itmap2) {
				pfl_name = itmap1->first;
				//pfl.setSens_fromShiftVals(shift_info_mp[0], itmap1->second, shift);
				pfl.setSens_fromShiftVals(shift_info_mp[1], itmap2->second, shift);
				//++itmap2;
			}
			pfl.print_sensitivities(dir_output + "sensitivities_total.csv");
			var_MC = pfl.calc_MC_Var(percentile);
					
			ff_var << "Top Portfolio," << var_MC << std::endl;
			pfl.print_vars(dir_var + "var_ir_report.csv");
			ff_var.close();	
		}
	}
	if(do_CalcVar) {
		if(!do_GenSims) { // need to read simulations from the file

		}
		percentile = params.get_percentile();
		var_MC = pfl.calc_MC_Var(percentile);
		pfl.print_vars(dir_var + "var_report.csv");
	}
	if(do_CCAR) {
		pfl.save_rets_eod();

		// Read CCAR info
		ccar_info.read_info(dir_CCAR);
		n_proj = ccar_info.getNumProjections();
		ccar_info.calc_crf_shifts(pfl, ccar_shifts);
		//pfl.digest_crf_info(ccar_info);
		
		for(ii = 0; ii < n_proj; ++ii) {
			pfl.apply_CCAR_info(ccar_shifts, ii, ccar_info.getFactTypeCur_coreFact_mp(), ccar_info.getReg_2_cur_mp(), cur_2_reg_map);
			
			std::ostringstream ostr;
			ostr << ii;
			pfl.output_sim_diff(dir_output + "ccar_sim_proj_" + ostr.str() + ".csv");

			//var_MC = pfl.calc_MC_Var(percentile);
			//pfl.print_vars(dir_var + "var_report.csv");
			pfl.clear_CCAR_info(ccar_info, ii);
			pfl.reset_rets_eod();
		}
		
	}
	return 0;
}
void imposeCorr_REFI_TURNOVER(std::vector<boost::shared_ptr<Factor> >& facs, double corr) {
	
	if(facs.size() < 2)
		return;
	std::vector<double> sig1(facs[0]->getVol()), sig2(facs[1]->getVol());
	matrix<double> r1(facs[0]->_sim_rets), r2(facs[1]->_sim_rets);
	if(sig1.size() > 0 && sig2.size() > 0 && sig1[0] > 0)
		r1 *= corr * sig2[0] / sig1[0];
	r2 *= sqrt(1.0 - corr * corr);
	facs[1]->_sim_rets = r1 + r2;
}
string set_MBS_Ad_Hoc_RF_params(string dir_hist_data, std::vector<boost::shared_ptr<Factor> >& facs, double& cor) {
	
	ifstream ff;
	ff.open(dir_hist_data + "MBS_Ad_Hoc_RF.csv");
	if(!ff.is_open())
		return std::string("File " + dir_hist_data + "MBS_Ad_Hoc_RF.cs is not found");

	string buf("");
	std::vector<std::string> values;
	double refi_std(0.106), turn_std(0.084), refi_eod(1.0), turn_eod(1.0), std(0.0), eod(1.0), one_day_scaler(sqrt(1.0 / 252.0)), scale;
	unsigned int i;
	std::vector<double> vec(1);
	cor = -0.5;
	getline(ff, buf); // Name, Mode, Stdev, Starting value, Corr
	for(i = 0; i < 2; ++i) {
		getline(ff, buf); 
		boost::algorithm::split(values, buf, boost::is_any_of(",")); 
		if(values.size() >= 5) {
			if(values[0] == "REFI_DAILY") {
				refi_std = atof(values[2].c_str());
				refi_eod = atof(values[3].c_str());
				cor = atof(values[4].c_str());
			}
			if(values[0] == "TURNOVER_DAILY") {
				turn_std = atof(values[2].c_str());
				turn_eod = atof(values[3].c_str());
				cor = atof(values[4].c_str());
			}
		}
	}
	for(i = 0; i < facs.size(); ++i) {
		if(facs[i]->getHistName() == "REFI_DAILY") {
			std = refi_std;
			eod = refi_eod;
		}
		if(facs[i]->getHistName() == "TURNOVER_DAILY") {
			std = turn_std;
			eod = turn_eod;
		}
		//vec[0] = std * one_day_scaler;
		vec = facs[i]->getVol();
		scale = std / vec[0];
		vec[0] = std;
		facs[i]->_sim_rets = facs[i]->_sim_rets * scale;
		facs[i]->setVol(vec);

		vec[0] = eod;
		facs[i]->setEODVal(vec);
	}
	return "";
}
double calc_Hist_Var(std::vector<boost::shared_ptr<Factor> >& factors, std::vector<date>& dts, double percentile) {
	std::vector<boost::shared_ptr<Factor> >::iterator it;
	unsigned int idate, ndates(dts.size()), m;
	double var(0.0);
	std::vector<double> pfl_pnl(ndates), fact_pnl(ndates);
	bool isPnLCalced(false);
	for(it = factors.begin(); it != factors.end(); ++it) {
		isPnLCalced = (*it)->calcPnL_HR(fact_pnl);
		if(isPnLCalced)
			for(idate = 0; idate < dts.size(); ++idate) {
				pfl_pnl[idate] += fact_pnl[idate];
			}
	}
	m = ndates * percentile - 1;
	sort(pfl_pnl.begin(), pfl_pnl.end());
	var = pfl_pnl[m];
	return var;
}

double genSimsOld(std::vector<boost::shared_ptr<Factor> >& f_notMapped, std::vector<boost::shared_ptr<Factor> >& f_Regressed, std::vector<boost::shared_ptr<Factor> >& f_Idio, 
	std::vector<boost::shared_ptr<Factor> >& f_PCA, std::vector<date>& dts_common, unsigned int nsim, std::vector<double>& wts, unsigned int sim_horizon, unsigned int randomness_control) {
	
	int errCode(-3000); // ndates = 0
	unsigned int ndates(dts_common.size()), i, nterms, j, k, ncol, nRets(0), n_idio(0), n_fact, theTermForMapping;
	boost::numeric::ublas::vector<double> rand(nsim), vec, shift, tilt;
	matrix<double> systematic_part;
	double scaler(1.0), vv(0), resid_vol;
	std::vector<unsigned int> terms, mapped_terms;
	matrix<double> randoms(nRets, nRets), all_scen;
	std::vector<boost::shared_ptr<Factor> >::iterator itfac, jtfac;
	
	n_fact = f_Regressed.size() * 2 + f_Idio.size() * 2;
	if(ndates > sim_horizon) 
		nRets = ndates - sim_horizon;
	n_fact += nRets;

	// Setup random numbers generator
	boost::mt19937                     gener(1);
    boost::normal_distribution<double> normal_dist;   // Normal Distribution
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal_dist);
	rng.engine().seed(1);
	//rng.distribution().reset(); // to change the seed

	// All factors without gaps are handled by nRets independent normal variables
	// Each Mapped and Idiosyncratic is handled by 2 independent normal variables to impose certain correlations across terms
	
	// Control randomness
	for(i = 0; i < randomness_control; ++i)
		vv = rng();

	randoms.resize(nsim, n_fact);
	for(k = 0; k < nsim; ++k) 
			for(i = 0; i < n_fact; ++i) 
				randoms(k, i) = rng();

	if(0)   // Standard deviation for each column is 1 after that
		gramm_shmidt_orth(randoms); // Orthogonalize columns

	// Generate simulated returns for the factors treated as idiosyncratic
	// If a factor has multiple terms it is simulated with 2 uncorrelated variables (PCA with 2 components sort of one shift and one tilt).
	ncol = nRets;

	for(itfac = f_Idio.begin(); itfac != f_Idio.end(); ++itfac) {
		std::vector<double> vols((*itfac)->getVol());
		nterms = vols.size();
		
		// Factors treated as idiosyncratic set with a single term
		(*itfac)->_sim_rets.resize(nsim, nterms);
		(*itfac)->calcAdhocShiftTilt(shift, tilt);

		for(i = 0; i < (*itfac)->getNterms(); ++i) 
			column((*itfac)->_sim_rets, i) = (column(randoms, ncol) * shift[i] + column(randoms, ncol + 1) * tilt[i])* vols[i]; 
		
		ncol += 2;
	}
	if(nRets == 0)
		return errCode;

	// Make variances equal to 1 so that randoms^T * randoms / nsim = E
	for(i = 0; i < n_fact; ++i) {
		scaler = inner_product(column(randoms, i).begin(), column(randoms, i).end(), column(randoms, i).begin(), 0.0) / nsim;
		scaler = 1.0 /sqrt(scaler);
		column(randoms, i) *= scaler;
	}

	// Simulated returns for the factors with "good" data are weighted combinations of rows of historical returns: R* = M * W^(-0.5) * R
	matrix_range<matrix<double> > mr(randoms, boost::numeric::ublas::range(0, nsim), boost::numeric::ublas::range(0, nRets));
	
	for(itfac = f_notMapped.begin(); itfac != f_notMapped.end(); ++itfac) 
		(*itfac)->genSims(mr, wts);
	
	for(itfac = f_PCA.begin(); itfac != f_PCA.end(); ++itfac) 
		(*itfac)->genSims(mr, wts);
	

	// Generate simulated returns for successfully mapped factors
	// A single term was regressed. For all others terms simulations are scaled to match historical volatility
	ncol = nRets + 2 * f_Idio.size();
	
	for(itfac = f_Regressed.begin(); itfac != f_Regressed.end(); ++itfac) {
		std::vector<double> vols((*itfac)->getVol());
		nterms = vols.size();
		(*itfac)->_sim_rets.resize(nsim, nterms);
		theTermForMapping = (*itfac)->getTheTermForMapping();

		mapped_terms = (*itfac)->_mapped_terms;
		matrix<double> sim_rets(nsim, mapped_terms.size());
		
		// Columns of matrix sim_rets are filled with simulated values for the "mapped to" factors
		j = 0;
		for(jtfac = (*itfac)->getMappedTo().begin(); jtfac != (*itfac)->getMappedTo().end(); ++jtfac) {
			column(sim_rets, j) = column((*jtfac)->_sim_rets, mapped_terms[j]);
			++j;
		}
		// Multiply by regression coef to get the matrix (actually the vector since only single term was regressed) and add idiosyncratic risk
		systematic_part = prod(sim_rets, (*itfac)->getRegrCoef());
		vec = column(systematic_part, 0);
		scaler = inner_product(vec.begin(), vec.end(), vec.begin(),0.0)/ nsim;
		scaler = 1.0 / sqrt(scaler);
		vec *= scaler;

		(*itfac)->calcAdhocShiftTilt(shift, tilt);
		resid_vol = (*itfac)->getResidVol();
		
		for(j = 0; j < nterms; ++j) 
			column((*itfac)->_sim_rets, j) = vols[j] * (vec + resid_vol * (column(randoms, ncol) * shift[j] +  column(randoms, ncol + 1) * tilt[j]));

		ncol += 2;;
	}
	return vv;
}
double calc_pca_projections(std::vector<boost::shared_ptr<Factor> >& facs, unsigned int n_eig, double decay, std::vector<boost::shared_ptr<Factor> >& f_PCA, 
	std::vector<boost::shared_ptr<Factor> >& f_PCA_all, unsigned int sim_horizon) {
		
	if(facs.size() == 0) // Should return error Tanya
		return 0.0;

	
	FactType factType(facs[0]->getType());
	std::vector<date> dts, dts_res, dts_common(facs[0]->getHrDates()), new_dts;
	if(dts_common.size() <= sim_horizon)
		return 0.0; // Should return error Tanya

	double sqrt_days_inv(1.0/sqrt(dts_common.size() * 1.0));
	date_duration dt1, dt2;

	// Calculates PCA for a set of factors. Returns the percentage of variance explained by these factots
	std::vector<double> wts;
	int errCode(0);
	unsigned int n, ndates(0), nrets(0), nfacs(facs.size()), nterms, i, j, maxDates(3600);
	double vv(0.0), acc, total_variance, percentage_explained(0.0);;
	matrix<double> eigvec, rets, rr, hr, hr_values;
	std::vector<double> eigval, vec(1);
	std::vector<unsigned int> terms(1, 0);
	matrix<double> sig_inv, eig_inv_sqr(n_eig, n_eig), m1, m2;
	
	ndates = dts.size();
	nrets = ndates - 1;

	// Find common dates
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::vector<date>::iterator it_res;
	for(itfac = facs.begin(); itfac != facs.end(); ++itfac) {
		dts = (*itfac)->getHrDates();
		dts_res.resize(maxDates);
		it_res = set_intersection(dts_common.begin(), dts_common.end(), dts.begin(), dts.end(), dts_res.begin());
		dts_res.resize(it_res - dts_res.begin());
		dts_common.assign(dts_res.begin(), dts_res.end());
	}
	ndates = dts_common.size();
	if(ndates < 2)
		return 0;

	nrets = ndates - sim_horizon;
	calculate_weights(decay, dts_common, wts, sim_horizon);

	// Collect returns for all factors
	for(itfac = facs.begin(); itfac != facs.end(); ++itfac) {
		(*itfac)->setRetsOnDates(dts_common, sim_horizon, (*itfac)->_rets, new_dts);
		rr = (*itfac)->_rets;
		nterms = rr.size2();
		n = rets.size2();
		rets.resize(nrets, n + nterms);
		
		for(j = 0; j < nterms; ++j)
			for(i = 0; i < nrets; ++i)
				rets(i, n + j) = rr(i, j) * wts[i];

	}
	nfacs = rets.size2();
	sig_inv.resize(nfacs, nfacs);
	sig_inv.clear();
	for(i = 0; i < nfacs; ++i) {
		vv = inner_product(column(rets, i).begin(), column(rets, i).end(), column(rets, i).begin(), 0.0);
		sig_inv(i, i) = 1.0 / sqrt(vv);
	}
	eig_info << factTypeNames[factType] << "," << nfacs << ",";
	if(nfacs >= n_eig) {
		svdLanczos(rets, n_eig, eigvec, eigval);
		n_eig = eigval.size();

		// Calculate historical returns (projections returns on principal components) to regress against: H = W^0.5 * R * sig^(- 0.5)  * Q  / sqrt(eigval), R has been already multiplied by W^0.5
		eig_inv_sqr.clear();
		acc = 0.0;
		for(i = 0; i < n_eig; ++i) {
			eig_inv_sqr(i,i) = 1.0/sqrt(eigval[i]);
			acc += eigval[i];
		}
		m1 = prod(eigvec, eig_inv_sqr);
		m2 = prod(sig_inv, m1);
		hr = prod(rets, m2);

		total_variance = nfacs;
		percentage_explained = acc / total_variance;
		
		for(i = 0; i < n_eig; ++i)
			eig_info << eigval[i] << ",";
		eig_info << percentage_explained << ",";
	}
	else {
		hr = prod(rets, sig_inv);
		n_eig = nfacs;
	}
	eig_info << n_eig << std::endl;

	double sim_horizon_days(1.0);
	std::vector<double> vals(dts_common.size());

	for(i = 0; i < n_eig; ++i) {
		boost::shared_ptr<Factor> ptr(new Factor);
		(*ptr).setType(facs[0]->getType());
		std::ostringstream ostr;
		ostr << i;
		if(factType == Eqt_Spot) {// To be consistent with existing naming conventions for equities
			if(i < 10)
				(*ptr).setHistName("J0" + ostr.str() + "IRVol");
			else
				(*ptr).setHistName("J" + ostr.str() + "IRVol");
			(*ptr).setScenName(":factorGB_" + ostr.str());
		}
		else {
			(*ptr).setHistName("PCA_" + factTypeNames[ptr->getType()] + string("_") + ostr.str());
			(*ptr).setScenName(":" + (*ptr).getHistName());
		}
		(*ptr).setNterms(1);
		(*ptr).setCurr("USD");
		(*ptr).setRegion("North America");
		(*ptr).setHrDates(dts_common);
		(*ptr).setR2(0.0);
		(*ptr).setResidVol(0.0);
		(*ptr).setTerms(terms);
		(*ptr).setMappingType(PCA);
		std::vector<double> eod(1), vol(1);
		std::vector<unsigned int> terms(1);
		terms[0] = 0;
		vol[0] = sqrt_days_inv; 
		eod[0] = 1.0;
		(*ptr).setVol(vol);
		(*ptr).setEODVal(eod);
		(*ptr).setTerms(terms);
		(*ptr).setProxyName("");
		for(j = 0; j < ndates; ++j) {
			if(j == 0)
				vals[j] = 1.0;
			else {
				dt1 = dts_common[sim_horizon] - dts_common[0];
				dt2 = dts_common[j] - dts_common[0];
				sim_horizon_days = dt1.days();
				if(j < sim_horizon)
					vals[j] = exp(hr(0, i) * dt2.days() / sim_horizon_days);
				else
					vals[j] = vals[j - sim_horizon] * exp(hr(j - sim_horizon, i));
			}
			vec[0] = vals[j];
			(*ptr)._hrDataMap.insert(std::pair<date, std::vector<double> > (dts_common[j], vec));
		}

		// scale to have variance of 1/sqrt(days)
		(*ptr)._rets.resize(nrets,1);
		for(j = 0; j < wts.size(); ++j) 
			(*ptr)._rets(j, 0) = sqrt_days_inv * hr(j, i) / wts[j];
	
		f_PCA.push_back(std::move(ptr));
		f_PCA_all.push_back(f_PCA[f_PCA.size() -1]);
	}
	/* // Write to file
		for(j = 0; j < eigval.size(); ++j) 
			fout << eigval[j] << ",";
		fout << endl;
		for(j = 0; j < hr.size1(); ++j) {
			for(i=0; i< hr.size2(); i++) {
				fout << hr(j,i) << ",";
			}
			fout << endl;
		}
		fout << endl;

		for(j = 0; j < eigvec.size1(); ++j) {
			for(i=0; i< eigvec.size2(); i++) {
				fout << eigvec(j,i) << ",";
			}
			fout << endl;
		}
	}
	fout.close();
	*/
	return percentage_explained;
}

int readHistData(ErrorLog& errLog, ParamControl& params, std::map<string, string>& cur_2_reg_map, string dir_hist_data, std::vector<string>& h_names, std::vector<boost::shared_ptr<Factor> >& factors) {
// Read given set of files with historical data and create array of factors
	
	string fn[]={"hr_curve.csv", "hr_spot.csv", "hr_spot_20140905_All.csv"}, dir_name; 
	//string fn[]={"hr_spot_20140526_small.csv"}; 
	static const std::vector<string> hr_file_names(fn, fn + sizeof(fn)/sizeof(fn[0]));
	string buf(""), fullName, thisName, factorGroup, curr;
	std::vector<string> values, factNames, dayMonYr;
	std::vector<date> hr_dates;
	std::vector<string>::iterator itvec;
	unsigned int numHrFiles(hr_file_names.size()), iFactor(0), icount(0), iterm(0), i, idate(0), day(1), year(1999);
	unsigned short month(1);
	int ln(0);
	std::ifstream inp_file;
	using boost::is_any_of;
	bool newFactorFound(false);
	double val;
	std::map<date, std::vector<double> > histDataMap;
	std::map<date, std::vector<unsigned int> >	 histTermsMap;
	std::vector<double> prc;
	std::vector<unsigned int> terms;
	std::vector<unsigned int> term0(1,0);
	date dd;
	std::ofstream file_trace;

	for(i=0; i < numHrFiles; ++i) { 
	//for(i=2; i < 3; ++i) { 
		// Todo: Comment out
		dir_name = dir_hist_data + hr_file_names[i];
		//dir_name =dir_hist_data + "clean_" + hr_file_names[i];
	
		inp_file.open(dir_name);
		if(!inp_file.is_open())
			continue;

		iFactor = 0;
		terms.clear();
		prc.clear();
		fullName = string("No Risk F with such name");
		hr_dates.clear();
		histDataMap.clear();
		histTermsMap.clear();
		
		// Todo: Comment out
		std::ofstream f_clean;
		f_clean.open(dir_hist_data + "clean_" +  hr_file_names[i]);
		std::vector<string> tmp;
		
		if(hr_file_names[i].find("hr_curve.csv") != string::npos || hr_file_names[i].find("hr_eqt_vol_curve.csv") != string::npos) { // reading curves with multiple terms
			
			while(getline(inp_file, buf)) {
				// Todo: Comment out
				tmp.push_back(buf.c_str());

				if (buf.length() == 0 ) 
					continue;
				// Can encounter the following strings:
				// 1) FactorGroup,Name,Currency,CurveType,CompoundingPeriod,DayCountBasis,Price/Yield,Date,Term,Value
				// 2) EQUITY_DATA,ACS_C00819010_EQTV,USD,,NA,NA,,2012/03/23,184,0.400152
				// 3) ,,,,,,,2012/03/26,93,0.2245
				// 4) ,,,,,,,,184,0.2311
				
				if(buf.find("FactorGroup") != string::npos) {  // new factor
					newFactorFound = true;
					idate = 0;
					continue;
				}

				boost::algorithm::split(values, buf, is_any_of(",")); 
				if(newFactorFound) { 

					itvec = std::find(h_names.begin(), h_names.end(), fullName);
					if(itvec != h_names.end()) { // Setup previous factor
						// save the last date
						if(prc.size() == terms.size() && prc.size() > 0) {
							histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
							histTermsMap.insert(std::pair<date, std::vector<unsigned int> >(dd, terms));
						}
						boost::shared_ptr<Factor> ptr(new Factor);
						ptr->setAttributes(cur_2_reg_map, "hr_curve", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
						ptr->dataCleanUp(params);
						factors.push_back(std::move(ptr));
						
						// Todo: Comment out
						
						ln = tmp.size();
						if(iFactor > 0)
							ln -= 2;
						for(unsigned int k = 0; k < ln; ++k)
							f_clean << tmp[k] << std::endl;
						
					}
					++iFactor;
					
					// Todo: Comment out
					if(tmp.size() > 1) 
						tmp.erase(tmp.begin(), tmp.end() - 2);
					
					newFactorFound = false;
					factorGroup = values[0].c_str();
					fullName = values[1].c_str();
					curr = values[2].c_str();

					if(fullName == "USDPrime" || fullName == "EURMM") 
						int tt = 0;

					histDataMap.clear();
					histTermsMap.clear();
					hr_dates.clear();
				}
				if(values.size() == 10) { 
					if(buf.find(",,,,,,,,") == string::npos ) { // new date
						// if it is not a new factor, save the previous date
						if(prc.size() == terms.size() && prc.size() > 0 && idate != 0) {
							histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
							histTermsMap.insert(std::pair<date, std::vector<unsigned int> >(dd, terms));
						}
						boost::algorithm::split(dayMonYr, values[7], is_any_of("/")); 
						year =atoi(dayMonYr[0].c_str());
						month =atoi(dayMonYr[1].c_str());
						day =atoi(dayMonYr[2].c_str());
						dd = date(year, greg_month(month), day);
						hr_dates.push_back(dd);

						prc.clear();
						terms.clear();
						++idate;
					}

					val = atof(values[9].c_str());
					int itrm = atoi(values[8].c_str());
					//if(val <= 0 || itrm <= 0) {
					if(itrm <= 0) {
						ostringstream c1, c2;
						c1 << val;
						c2 << itrm;
						string line(fullName + ", not valid term," + values[7] + "," + c1.str() + "," + c2.str());
						errLog.add_line(line);
					}
					else {
						prc.push_back(val);
						terms.push_back(itrm);
					}
				}
			} // end while
			// Last factor
			itvec = std::find(h_names.begin(), h_names.end(), fullName);
			if(itvec != h_names.end()) { 
				boost::shared_ptr<Factor> ptr(new Factor);
				
				if(prc.size() == terms.size() && prc.size() > 0) {
					histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
					histTermsMap.insert(std::pair<date, std::vector<unsigned int > >(dd, terms));
				} 
				ptr->setAttributes(cur_2_reg_map, "hr_curve", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
				ptr->dataCleanUp(params);
				factors.push_back(std::move(ptr));
				
				// Todo: Comment out
				for(unsigned int k = 0; k < tmp.size(); ++k)
					f_clean << tmp[k] << std::endl;
				
			}
			// Todo: Comment out
			tmp.clear();
		}
		else { // reading single term data
			
			while(getline(inp_file, buf)) {
				// Todo: Comment out
				tmp.push_back(buf.c_str());
				
				if (buf.length() == 0)
					continue;
			    
				if(buf.find("FactorGroup") != string::npos) {  // new factor
					newFactorFound = true;
					idate = 0;
					continue;
				}

				boost::algorithm::split(values, buf, is_any_of(",")); 
				if(newFactorFound) {
					itvec = std::find(h_names.begin(), h_names.end(), fullName);
					if(itvec != h_names.end()) { // set up previous factor
						boost::shared_ptr<Factor> ptr(new Factor); 
						ptr->setAttributes(cur_2_reg_map, "hr_spot", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
						ptr->dataCleanUp(params);
						factors.push_back(std::move(ptr));
						
						// Todo: Comment out
						
						ln = tmp.size();
						if(iFactor > 0)
							ln -= 2;
						for(unsigned int k = 0; k < ln; ++k)
							f_clean << tmp[k] << std::endl;
						
					}
					++iFactor;
					
					// Todo: Comment out
					if(tmp.size() > 1) 
						tmp.erase(tmp.begin(), tmp.end() - 2);
					
					newFactorFound = false;
					factorGroup = values[0].c_str();
					fullName = values[1].c_str();
					if(fullName == "C00036020_NASDAQ_AAONSpot" || fullName == "C00036110_NYSE_AIRSpot") {
						int tt = 0;
					}
					curr = values[2].c_str();
					hr_dates.clear();
					histDataMap.clear();
					histTermsMap.clear();
				}
				
				if(values.size() == 6) {
					boost::algorithm::split(dayMonYr, values[4], is_any_of("/")); 
					year =atoi(dayMonYr[0].c_str());
					month =atoi(dayMonYr[1].c_str());
					day =atoi(dayMonYr[2].c_str());
					dd = date(year, greg_month(month), day);

					val = atof(values[5].c_str());
						hr_dates.push_back(dd);
						prc.push_back(val);
						histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
						histTermsMap.insert(std::pair<date, std::vector<unsigned int> >(dd, term0));
						prc.clear();

					++idate;
				}
			} // end while
			itvec = std::find(h_names.begin(), h_names.end(), fullName);
			if(itvec != h_names.end()) {
				boost::shared_ptr<Factor> ptr(new Factor); 
				ptr->setAttributes(cur_2_reg_map, "hr_spot", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
				ptr->dataCleanUp(params);
				factors.push_back(std::move(ptr));
				
				// Todo: Comment out
				for(unsigned int k = 0; k < tmp.size(); ++k)
					f_clean << tmp[k] << std::endl;
			}
			// Todo: Comment out
			tmp.clear();
		}
		inp_file.close();
		// Todo: Comment out
		f_clean.close();
	}

	return 0;
}

int readScenarioFile(string inp_fn, std::vector<boost::shared_ptr<Factor> >& factors, unsigned int nsim) {
	unsigned int maxTermsNum(100), idx, isim;
	std::ifstream inp_file;
	inp_file.open(inp_fn);
	std::ofstream file_factors; 
	string fact_info_fn("G:\\VAR attribution\\20140224\\Output\\fact_info.csv");
	file_factors.open(fact_info_fn);
	//matrix<double> sim_rets(nsim, maxTermsNum);
	Factor*  factor_cur;
	std::vector<unsigned int> terms;
	string buf("");
	std::vector<string> values;
	std::vector<string> nameParts;
	std::vector<double>  prc(maxTermsNum, 0);  
	string curName;
	int unIdentified(0), ic1(0), ic2(0), ic3(0), ic4(0), ic6(0);
	int scenNum(-1), iFactor(0);
	
	// Read scenario file to get the set of factors
	while(getline(inp_file, buf)) {
		if (buf.length() == 0)
			continue;
		if(buf.find("base")!=string::npos) 
				break;
		if (buf.find("Scenario Set") == string::npos){ // start new scenario
			if (buf.find("Mmc2000")!=string::npos) {
				++scenNum;
				iFactor = 0;
				if(scenNum == 3)
					break;
			}
			
			if (buf.find("msmc 2000")!=string::npos || buf.find("green")!=string::npos ||
					buf.find("Constant")!=string::npos || buf.find("Trigger")!=string::npos ||
					buf.find("Term")!=string::npos) { // encountered new factor
				
				if((scenNum == 0 && iFactor > 0) || (scenNum ==1 && iFactor == 0)) {
					boost::shared_ptr<Factor> ptr(factor_cur);
					ptr->setNterms(terms.size());
					ptr->setTerms(terms);
					ptr->setScenName(curName);
					string h_name = s_name_2_h_name(curName);
					ptr->setHistName(h_name);

					if(h_name == "USDInterbank")
						int tt = 0;

					ptr->factPrint(file_factors);
					ptr->setMappingType(NotMapped);
					ptr->_sim_rets.resize(nsim, terms.size());
					factors.push_back(std::move(ptr));
				}
				if( !(iFactor == 0 && scenNum == 0)) {
					idx = iFactor;
					isim = scenNum;
					if(iFactor == 0 && scenNum > 0) {
						idx = factors.size();
						isim = scenNum -1;
					}
					for(unsigned j = 0; j < terms.size(); ++j) 
						factors[idx - 1]->_sim_rets(isim, j) = prc[j];
				}
				terms.clear();
				
				using boost::is_any_of;                                                                                                                                     
				boost::algorithm::split(values, buf, is_any_of(",")); 
				curName = values[5];
				//if(scenNum > 0)
				//	break;
				
				boost::algorithm::split(nameParts, curName, is_any_of("_")); 
				factor_cur = new Factor;
				populateFactor(factor_cur, nameParts, ic1, ic2, ic3, ic4, ic6, unIdentified);
				++iFactor;
			}
			else if(values.size() == 15) { // read terms
				boost::algorithm::split(values, buf, is_any_of(",")); 
				prc[terms.size()] = atof(values[14].c_str());
				terms.push_back(atoi(values[13].c_str()));	
			}
		}

	} // end while

	inp_file.close();
	file_factors.close();

	return 0;
}

int alignDates(ErrorLog& errLog, ParamControl& params, std::vector<boost::shared_ptr<Factor> >& factors, std::vector<date>& dts_common) {
	unsigned int maxDates(params.get_maxDates()), minDatesCount(params.get_minDatesCount()), minNumRetsForRegr(params.get_minNumRetsForRegr()), i, j, sim_horizon(params.getSimHorizon());
	std::vector<date> dts, dts_res;
	double defaultIdioVol(params.get_defaultIdioVol()), defaultIdioEquityVol(params.get_defaultIdioEquityVol()), decay(params.get_decay_factor());
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::vector<date>::iterator it_res;
	dts_res.resize(maxDates);
	std::map<date, std::vector<double> >::iterator it;
	std::vector<string>::const_iterator itstr;
	
	if(factors.size() == 0)
		return 0;

	// Sort factors that have sufficiently long history and positive prices from the rest and identify the set of dates common across these factors
	i = 0;
	itfac = factors.begin();
	while(itfac != factors.end() &&  (*itfac)->getHrDates().size() < minDatesCount) {
		(*itfac)->setDataAvalType(WithGaps);
		++itfac;
	}
	if(itfac == factors.end())
		--itfac;

	dts_common = (*itfac)->getHrDates();
	while(itfac != factors.end()) {

		if((*itfac)->getHistName() == "C00036020_NASDAQ_AAONSpot")
			int tt = 0; // Temp Tanya

		dts = (*itfac)->getHrDates();
		dts_res.resize(maxDates);
		it_res = set_intersection(dts_common.begin(), dts_common.end(), dts.begin(), dts.end(), dts_res.begin());
		dts_res.resize(it_res - dts_res.begin());
		
		if(dts.size() >= dts_common.size() && dts_res.size() < minDatesCount )  {// data completion is needed
			(*itfac)->complete(dts_common);
			dts_res = dts_common;
		}
		if(dts_res.size() < minDatesCount) {  
			(*itfac)->setDataAvalType(WithGaps);
		}
		else	{
			(*itfac)->setDataAvalType(NoGaps);
			dts_common.assign(dts_res.begin(), dts_res.end());
			j = dts_common.size();
			
		}
		++itfac;
	}

	// Setup historical data, returns and vols on the common dates
	for(itfac=factors.begin(); itfac != factors.end(); ++itfac) 
		if((*itfac)->getDataAvalType() == NoGaps) {
			(*itfac)->setMappingType(NotMapped);
			(*itfac)->setRetsVols(params, dts_common);
		}

	return 0;
}

int alignSens_HistFactors(std::vector<boost::shared_ptr<Factor> >& s_factors, std::vector<boost::shared_ptr<Factor> >& h_factors) {
	// Assign spot and historical data to the sensitivity factors
	std::vector<boost::shared_ptr<Factor> >::iterator it, jt;
	string factName;
	std::ofstream assign_fact_log;
	assign_fact_log.open("G:\\VAR attribution\\20140224\\Output\\log_assign.csv");

	for(it = s_factors.begin(); it != s_factors.end(); ++it) {
		factName = (*it)->getScenName();
		jt = h_factors.begin(); 
		
		while(jt != h_factors.end() && (*it)->comp_Scen_Hist_names((*jt)->getCurr(), (*jt)->getSubType(), (*jt)->getType(), (*jt)->getHistName()) != 0)
			++jt;
		if(jt != h_factors.end()) {
			//(*it)->setRets((*jt)->getRets().get());
			(*it)->_hrDataMap = (*jt)->_hrDataMap;
			(*it)->setHrDates((*jt)->getHrDates());
			if((*it)->getTerms().size() != (*jt)->getTerms().size())
				assign_fact_log << "alignSens_HistFactors: Different terms for historical and sensitivity data: " << factName << endl;
			//(*it)->setSpot((*jt)->getSpot());
			//(*jt)->resetRets();
			
		}
		else
			assign_fact_log << "alignSens_HistFactors: No hist data found for: " << factName << endl;
	}
	assign_fact_log.close();
	return 0;
}

unsigned int getNumFactors(std::vector<boost::shared_ptr<Factor> >& f_vec) {
	unsigned int n(0);
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = f_vec.begin(); itfac != f_vec.end(); ++itfac) 
		n += (*itfac)->getNterms();
	
	return n;
}

int read_cur_2_reg_map(string map_fn, std::map<string, string>& cur_2_reg_map) {
	string buf(""), factName, cur, subType, reg;
	std::vector<string> values;
	int errCode(-1001);

	std::ifstream cur_2_reg;
	cur_2_reg.open(map_fn);

	using boost::is_any_of; 
	getline(cur_2_reg, buf);
	boost::algorithm::split(values, buf, is_any_of(",")); 
	if(values.size() != 2)
		return errCode;

	while(getline(cur_2_reg, buf)) {
		boost::algorithm::split(values, buf, is_any_of(",")); 
		if(values.size() != 2)
			return errCode;
		cur = values[0].c_str();
		reg = values[1].c_str();

		cur_2_reg_map.insert(std::pair<string, string>(cur, reg));
	}

	cur_2_reg.close();
	return 0;
}
// This should be made as a template to invert any map
int invert_cur_2_reg(std::map<string, string>& cur_2_reg_map, std::map<string, std::vector<string> >& reg_2_cur_map) {
	std::map<string, string>::iterator itinp;
	std::map<string, std::vector<string> >::iterator itout;
	int errCode(0);
	string reg, cur;
	std::vector<string> currencies;

	for(itinp = cur_2_reg_map.begin(); itinp != cur_2_reg_map.end(); ++itinp) {
		cur = itinp->first;
		reg = itinp->second;
		itout = reg_2_cur_map.find(reg);
		
		currencies.clear();
		
		if(itout != reg_2_cur_map.end())	{
			currencies = itout->second;
			reg_2_cur_map.erase(itout);

		}
		currencies.push_back(cur);
		reg_2_cur_map.insert(std::pair<string, std::vector<string> >(reg, currencies));
	}

	return errCode;
}

void read_RF_names(string dir, std::vector<string>& h_names) {
	// The lists of risk factors names for yesterdays positions are under:
	// /prod/vol2/alg/batch/RWH/20140501/input/marketData

	//	For equity: The string "Spot" should be appended to the end of the names in this file in order to match the names in files with historic data
	// eqtRWHOpenPositionReport.csv
	
	string buf(""), vv, s_name, h_name;
	std::ifstream inp_file;

	// Use unique list of names prepared previously
	/*
	inp_file.open(dir + "unique_rf_names.csv");
	while(getline(inp_file, buf)) 
			h_names.push_back(buf.c_str());
		
	inp_file.close();
	return;
	*/
	// Create the unique set of names
	string f_list[] = {string("riskFactorTemplateIR.csv"), string("riskFactorTemplateCMSR.csv"), string("riskFactorTemplateEQUITY.csv"), string("riskFactorTemplateBNYE.csv"), 
			string("riskFactorTemplateBNYMCMNOSPRISK.csv"), string("riskFactorTemplateBNYMCMSPRISK.csv"), 
			string("riskFactorTemplateCDS.csv"), string("riskFactorTemplateCDSSR.csv"), string("riskFactorTemplateCMNOSR.csv"),  
			string("riskFactorTemplateCstIR.csv"), string("riskFactorTemplateEQUITYSR.csv"), string("riskFactorTemplateFX.csv"), 
			string("riskFactorTemplateFXNEW.csv"), string("riskFactorTemplatePERSHING.csv"), string("riskFactorTemplatePERSHINGNOSR.csv"), 
			string("riskFactorTemplatePERSHINGSR.csv")};
	//string eqt_list[] = {string("fromCmiMurexPrevDay.csv")}; 
	string eqt_list[] = {string("eqtRWHOpenPositionReport.csv"), string("cmiequityCurrIssues.csv")}; // to get 4856 factors 
	//string eqt_list[] = {string("fromFactorCoefGB_rerunable.csv")}; // to get 6909 stocks
	//string eqt_list[] = {string("fromGbFactorCoefOri.csv")};  // to get around 15000 stocks
	std::vector<string> fn(f_list, f_list + sizeof(f_list)/sizeof(f_list[0])), fn_eq(eqt_list, eqt_list + sizeof(eqt_list)/sizeof(eqt_list[0]));
	std::vector<string> values, names_vec, names;
	
	unsigned int n(fn.size()), i;
	std::vector<string>::iterator it;
	
	for(i = 0; i < n; ++i) { 
	//for(i = 0; i < 1; ++i) {
		inp_file.open(dir + fn[i]);
		while(getline(inp_file, buf)) {
			if(buf == ":factorGB_0" || buf == ":factorGB_1" || buf == ":factorGB_2" || buf == ":factorGB_3" || buf == ":factorGB_4" || buf == ":factorGB_5" || buf == ":factorGB_6" || buf == ":factorGB_7" || 
				buf == ":factorGB_8" || buf == ":factorGB_9" || buf == ":factorGB_10" || buf == ":factorGB_11" || buf == ":factorGB_12" || buf == ":factorGB_13" || buf == ":factorGB_14" || 
				buf == ":factorGB_15" || buf == ":factorGB_16" || buf == ":factorGB_17") // These are PCA and will be created
				continue;

 			vv = s_name_2_h_name(buf);
			names.push_back(vv);
		}
		inp_file.close();
	}
	
	for(i = 0; i < fn_eq.size(); ++i) {
		inp_file.open(dir + fn_eq[i]);
		while(getline(inp_file, buf)) {
			vv = buf + "Spot";
			// This is because of inconsistent naming convention Temp Tanya
			if(buf == "C85299600_Index_SPTR")
				vv = string("C85299600_OTC_SPTRSpot");
			if(buf == "C12496G10_Index_RUSSELL3000")
				vv = string("C12496G10_OTC_RUSSELL3000Spot");
			if(buf == "C12497K10_Index_VIX")
				vv = string("C12497K10_OTC_VIXSpot");
			if(buf == "C00000117_Index_SP500")
				vv = string("C00000117_Index_SPALNSpot"); // that's correct
			if(buf == "C26099405_Index_INDUA")
				vv = string("C26099405_OTC_INDUASpot");
			names.push_back(vv);
		}
		inp_file.close();
	}
	
// Make unique
	std::sort(names.begin(), names.end());
	it = std::unique(names.begin(), names.end());
	names.resize(std::distance(names.begin(), it));
	
	std::ofstream f_out;
	f_out.open(dir + "unique_rf.csv");
	for(i = 0; i < names.size(); ++i) {
		f_out << names[i] << std::endl;
		if(names[i].substr(0, 3) != "ZZZ")
			h_names.push_back(names[i]);
	}
	f_out.close();
	
	return;
}
 
 
int read_scen(std::string f_name, int nScen, std::vector<string>& names, matrix<double>&mat ) {
	// Reads scenarios that are in matrix form
	std::vector<string> values;
	string buf("");
	std::ifstream file;
	file.open(f_name);

	unsigned int icount(0), j(0);
	while(getline(file, buf)) {
		                                                                                                                                    
		boost::algorithm::split(values, buf, is_any_of(","));
		if(icount) {
			for(unsigned int i = 0; i < j; ++i)
				mat(icount - 1, i) = atof(values[i].c_str());
		}
		else {
			while(j < values.size()) {
				names.push_back(values[j]);
				++j;
			}
			mat.resize(nScen, j);
		}
		++icount;
	}

	return 0;
}
void assignCCAR_Mapping(std::vector<boost::shared_ptr<Factor> > &factors) {
	
	std::vector<boost::shared_ptr<Factor> >::iterator it;
	for(it = factors.begin(); it != factors.end(); ++it) {
		(*it)->setForecast();
	}

}
void populateFactor(Factor* factor, std::vector<string> nameParts, int& ic1, int& ic2, int& ic3, int& ic4, int& ic6, int& unIdentified) {

	FactType type(UnIdentified);
	double factVal(0.0);
	int term(0);
	string subType;
	string curr;
	bool identified(false);
	bool found(true);

	switch (nameParts.size()) {
		case 1: //FX spot: FXAED, 
			type = FX_Spot;
			curr = nameParts[0].substr(3,3);
			++ic1;
			break;
		case 2:  // IR curve: IREUR_Mmarket_304, gbFactors: factorGB_1_0, turnover: TURNOVER_DAILY_0, REFI: REFI_DAILY_0, CAP vol: index_cap3m11_0, Swaption vol: index_euro1_0,
			// FX vol: AUD_USD-Simulation_0, 
			if(nameParts[0] == ":factorGB") {
				type = factor_GB;
				curr = "USD";
			}
			else if(nameParts[0] == ":TURNOVER") {
					type = Turnover;
					curr = "USD";
				}
				else if(nameParts[0] == ":REFI") {
						type = Refi;
						curr = "USD";
					}
					else if(nameParts[0] == ":index") {
							if(nameParts[1].compare(0, 3, "cap")==0) {
								type = IRCap_Vol;
								curr = "USD";
							}
							else if(nameParts[1].compare(0, 4, "euro")==0) {
								type = IRSwaption_Vol;
								curr = "USD";
							}
							else
								found = false;
						}
					
						else if(nameParts[1] == "Forward") {
							type = FX_Forward;
							curr = nameParts[0].substr(3,3);
						}
						else if(nameParts[0].compare(0, 3, ":IR") == 0) {
								type = IR_Curve;
								curr = nameParts[0].substr(3,3);
								subType = nameParts[1];
							}
							else if(nameParts[1].size() > 10) {
								int len = nameParts[1].size();
								if(nameParts[1].compare(len-10, 10, "Simulation") == 0) {
									type = FX_Pair_Vol;
									curr = nameParts[0].substr(1,3) + "_"+nameParts[1].substr(0,3);
								}
							}
							
						if(found) 
							++ic2;
						else
							++unIdentified;

				break;
			case 3: // EQ_Vol: :AAPL_C03783310_EQTV
				if (nameParts[2] == "EQTV") {
					type = EQ_Vol;
					curr = "USD";
					++ic3;
				}
				else
					++unIdentified;
				break;
			case 4: // MBS: IRUSD_MBS_CC_OAS_0
				if(nameParts[1] == "MBS") {
					type = MBS;
					curr = "USD";
					++ic4;
				}
				else
					++unIdentified;
				break;
			case 5: // CDS: :CDS_AMERAXEL_UU2679_SNRFOR_EQTV
				if(nameParts[0] == ":CDS") {
					type = CDS;
					curr = "USD";
					++ic6;
				}
				else
					++unIdentified;
				break;
			case 6: // CDS: :CDS_ABBEYNATFP_GNF882_SNRFOR_MM_EQTV_365
				if(nameParts[0] == ":CDS") {
					type = CDS;
					curr = "USD";
					++ic6;
				}
				else
					++unIdentified;
				break;
			default:
				++unIdentified;
		}
		factor->setType(type);
		factor->setCurr(curr);
		factor->setSubType(subType);
		
		return;
}

string readShiftScen_Values(std::string fn, std::map<std::string, std::vector<double> >& pfl_val_map) {
	// Read valuation results from the file
	std::ifstream val_file;
	val_file.open(fn);
	string result;

	if(!val_file.is_open() )
		return string("readShiftScen_Values: File " + fn + " is not open");

	std::string buf(""), name;
	unsigned int nscen, nport, i, iscen;
	double val;
	std::vector<std::string> values;
	std::vector<double> prc;

	// Read from valuation file
	// First line: Number of Scenarios, Number of FactSets
	// Port1 0 value0
	// Port1 1 value1
	// ...
	// Port2 0 value0
	// Port2 1 value1
	getline(val_file, buf);
	boost::algorithm::split(values, buf, boost::is_any_of(",")); 
	if(values.size() != 3) 
		return "readShiftScen_Values: The number of fields in the header is not equal to 3";

	nscen = atoi(values[0].c_str());
	nport = atoi(values[1].c_str());

	for(i = 0; i < nport; ++i) {
		prc.clear();
		for(iscen = 0; iscen < nscen; ++iscen) {
			if(!getline(val_file, buf))
				return string("readShiftScen_Values: Empty line in " + fn);
			boost::algorithm::split(values, buf, boost::is_any_of(",")); 
			if(values.size() != 3) 
				return string("readShiftScen_Values: The number of fields in the line is not equal to 3, file " + fn);
			name = values[0].c_str();
			val = atof(values[2].c_str());
			prc.push_back(val);
		}
		pfl_val_map.insert(std::make_pair(name, prc));
	}
	return "";
}
string readShiftScen_FactorInfo(std::string info_fn, std::vector<std::string>& fact_names, std::vector<unsigned int>& nterms, double& shift) {
	// Read from the file the information about risk factors (RF) that was used to generate shift scenarios for sensitivities
	// This information is used to assigne sensitivities calculated from caluation files data to the RF
	std::ifstream info_file;
	info_file.open(info_fn);
	bool badresult(false);
	string result("");

	if(!info_file.is_open())
		return string("readShiftScen_FactorInfo: file " + info_fn + " can't be open");

	std::string buf(""), name_prev("newName"), name;
	unsigned int i(0), n;
	std::vector<std::string> values;
	
	// Read from factor info file
	getline(info_file, buf); // comment, shift, shift size
	boost::algorithm::split(values, buf, boost::is_any_of(",")); 
	if(values.size() != 3)
		return string("readShiftScen_FactorInfo: the number of fields in the header of " + info_fn + " is not equal to 3");

	shift = atof(values[2].c_str());

	getline(info_file, buf); // header ( "scenario,factor name, tenor")
	while(getline(info_file, buf)) {
		boost::algorithm::split(values, buf, boost::is_any_of(",")); 

		name = values[1].c_str();
		if(name_prev != name) {
			if(i != 0) {
				fact_names.push_back(name_prev);
				if(values.size() == 5)  // scenname, fact scen name, fact hist name, tenor, val
					n = i;
				else if(values.size() == 7) // Cross-Gammas scen name, fact scen name, fact hist name, tenor1, ternor2, val1, val2
					n = (1 + sqrt(1.0 + 8.0 * i)) * 0.5;
				else
					return string("readShiftScen_FactorInfo: the number of fields in the header of " + info_fn + " is not equal to 4 or 6");;
				nterms.push_back(n);
			}
			i = 0;
		}
		++i;
		name_prev = name;
	}
	if(name_prev != "newName") {
		fact_names.push_back(name_prev);
		nterms.push_back(n);
	}
	return result;
}

void calc_hr_eqt_factors_and_coef(std::vector<boost::shared_ptr<Factor> >& facs, double decay, unsigned int sim_horizon, std::string fn_hr, std::string fn_coef, std::vector<boost::shared_ptr<Factor> >& f_Mapped) {
	// Built the orthomormal basis in the historical returns space, calculate projections of all historical returns (for factors with and without missing data) to the basis
	if(facs.size() <1)
		return;

	unsigned int j(0), k, i, count_basis(0), nrets(facs[0]->_rets.size1()), nfacs(facs.size()), n(min(nrets, nfacs));
	std::vector<date> dts(facs[0]->getHrDates());
	std::vector<double> wts;
	matrix<double> basis(nrets, n), vec(nrets, 1), coef(n, 1);
	double sum(0), norm(0.0), min_norm(0.0001);
	std::vector<boost::shared_ptr<Factor> > f_Basis;

	calculate_weights(decay, dts, wts, sim_horizon);
	for(i = 0; i < nrets; ++i)
		wts[i] = wts[i] * wts[i];

	if(wts.size() != nrets)
		return;

	// Find basis
	while(j < nfacs && count_basis < n) {
		column(vec, 0) = column(facs[j]->_rets, 0);

		// Project next equity factor on the basis
		for(k = 0; k < count_basis; ++k) {
			sum = 0.0;
			for(i = 0; i < nrets; ++i)
				sum += wts[i] * basis(i, k) * vec(i, 0);
			coef(k, 0) = sum;
		}
		for(k = 0; k < count_basis; ++k) 
			column(vec, 0) -= coef(k, 0) * column(basis, k);
			
		norm = 0.0;
		for(i = 0; i < nrets; ++i)
			norm += vec(i, 0) *vec(i, 0) * wts[i];
		norm = sqrt(norm);
		if(norm > min_norm) {
			norm = 1.0 / norm;
			column(basis, j) = column(facs[j]->_rets, 0) * norm;
			++count_basis;
		}
		++j;
	}
	// Project on basis
	for(j = 0; j < facs.size(); ++j) {
		// Project next equity factor on the basis
		column(vec, 0) = column(facs[j]->_rets, 0);
		for(k = 0; k < count_basis; ++k) {
			sum = 0.0;
			for(i = 0; i < nrets; ++i)
				sum += wts[i] * basis(i, k) * vec(i, 0);
			coef(k, 0) = sum;
		}
		facs[j]->setRegrCoef(coef);
	}
	for(j = 0; j < f_Mapped.size(); ++j) {
		// This needs to be implemented
	}
}

int readHistData_clean(ErrorLog& errLog, ParamControl& params, std::map<string, string>& cur_2_reg_map, string dir_hist_data, std::vector<string>& h_names, std::vector<boost::shared_ptr<Factor> >& factors) {
// Read given set of files with historical data and create array of factors
	
	string fn[]={"hr_curve.csv", "hr_spot.csv", "hr_spot_20140905_All.csv"}, dir_name; 
	//string fn[]={"hr_spot_20140526_small.csv"}; 
	static const std::vector<string> hr_file_names(fn, fn + sizeof(fn)/sizeof(fn[0]));
	string buf(""), fullName, thisName, factorGroup, curr;
	std::vector<string> values, factNames, dayMonYr;
	std::vector<date> hr_dates;
	std::vector<string>::iterator itvec;
	unsigned int numHrFiles(hr_file_names.size()), iFactor(0), icount(0), iterm(0), i, idate(0), day(1), year(1999);
	unsigned short month(1);
	int ln(0);
	std::ifstream inp_file;
	using boost::is_any_of;
	bool newFactorFound(false);
	double val;
	std::map<date, std::vector<double> > histDataMap;
	std::map<date, std::vector<unsigned int> >	 histTermsMap;
	std::vector<double> prc;
	std::vector<unsigned int> terms;
	std::vector<unsigned int> term0(1,0);
	date dd;
	std::ofstream file_trace;

	for(i=0; i < numHrFiles; ++i) { 
	//for(i=0; i < 1; ++i) { 
		// Todo: Comment out
		//dir_name = dir_hist_data + hr_file_names[i];
		dir_name =dir_hist_data + "clean_" + hr_file_names[i];
	
		inp_file.open(dir_name);
		if(!inp_file.is_open())
			continue;

		iFactor = 0;
		terms.clear();
		prc.clear();
		fullName = string("No Risk F with such name");
		hr_dates.clear();
		histDataMap.clear();
		histTermsMap.clear();
		
		// Todo: Comment out
		//std::ofstream f_clean;
		//f_clean.open(dir_hist_data + "clean_" +  hr_file_names[i]);
		//std::vector<string> tmp;
		
		if(hr_file_names[i].find("hr_curve.csv") != string::npos || hr_file_names[i].find("hr_eqt_vol_curve.csv") != string::npos) { // reading curves with multiple terms
			
			while(getline(inp_file, buf)) {
				// Todo: Comment out
				//tmp.push_back(buf.c_str());

				if (buf.length() == 0 ) 
					continue;
				// Can encounter the following strings:
				// 1) FactorGroup,Name,Currency,CurveType,CompoundingPeriod,DayCountBasis,Price/Yield,Date,Term,Value
				// 2) EQUITY_DATA,ACS_C00819010_EQTV,USD,,NA,NA,,2012/03/23,184,0.400152
				// 3) ,,,,,,,2012/03/26,93,0.2245
				// 4) ,,,,,,,,184,0.2311
				
				if(buf.find("FactorGroup") != string::npos) {  // new factor
					newFactorFound = true;
					idate = 0;
					continue;
				}

				boost::algorithm::split(values, buf, is_any_of(",")); 
				if(newFactorFound) { 

					itvec = std::find(h_names.begin(), h_names.end(), fullName);
					if(itvec != h_names.end()) { // Setup previous factor
						// save the last date
						if(prc.size() == terms.size() && prc.size() > 0) {
							histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
							histTermsMap.insert(std::pair<date, std::vector<unsigned int> >(dd, terms));
						}
						boost::shared_ptr<Factor> ptr(new Factor);
						ptr->setAttributes(cur_2_reg_map, "hr_curve", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
						ptr->dataCleanUp(params);
						factors.push_back(std::move(ptr));
						
						// Todo: Comment out
						/*
						ln = tmp.size();
						if(iFactor > 0)
							ln -= 2;
						for(unsigned int k = 0; k < ln; ++k)
							f_clean << tmp[k] << std::endl;
						*/
					}
					++iFactor;
					
					// Todo: Comment out
					//if(tmp.size() > 1) 
					//	tmp.erase(tmp.begin(), tmp.end() - 2);
					
					newFactorFound = false;
					factorGroup = values[0].c_str();
					fullName = values[1].c_str();
					curr = values[2].c_str();

					if(fullName == "USDPrime" || fullName == "EURMM") 
						int tt = 0;

					histDataMap.clear();
					histTermsMap.clear();
					hr_dates.clear();
				}
				if(values.size() == 10) { 
					if(buf.find(",,,,,,,,") == string::npos ) { // new date
						// if it is not a new factor, save the previous date
						if(prc.size() == terms.size() && prc.size() > 0 && idate != 0) {
							histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
							histTermsMap.insert(std::pair<date, std::vector<unsigned int> >(dd, terms));
						}
						boost::algorithm::split(dayMonYr, values[7], is_any_of("/")); 
						year =atoi(dayMonYr[0].c_str());
						month =atoi(dayMonYr[1].c_str());
						day =atoi(dayMonYr[2].c_str());
						dd = date(year, greg_month(month), day);
						hr_dates.push_back(dd);

						prc.clear();
						terms.clear();
						++idate;
					}

					val = atof(values[9].c_str());
					int itrm = atoi(values[8].c_str());
					//if(val <= 0 || itrm <= 0) {
					if(itrm <= 0) {
						ostringstream c1, c2;
						c1 << val;
						c2 << itrm;
						string line(fullName + ", not valid term," + values[7] + "," + c1.str() + "," + c2.str());
						errLog.add_line(line);
					}
					else {
						prc.push_back(val);
						terms.push_back(itrm);
					}
				}
			} // end while
			// Last factor
			itvec = std::find(h_names.begin(), h_names.end(), fullName);
			if(itvec != h_names.end()) { 
				boost::shared_ptr<Factor> ptr(new Factor);
				
				if(prc.size() == terms.size() && prc.size() > 0) {
					histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
					histTermsMap.insert(std::pair<date, std::vector<unsigned int > >(dd, terms));
				} 
				ptr->setAttributes(cur_2_reg_map, "hr_curve", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
				ptr->dataCleanUp(params);
				factors.push_back(std::move(ptr));
				
				// Todo: Comment out
				//for(unsigned int k = 0; k < tmp.size(); ++k)
				//	f_clean << tmp[k] << std::endl;
				
			}
			// Todo: Comment out
			//tmp.clear();
		}
		else { // reading single term data
			
			while(getline(inp_file, buf)) {
				// Todo: Comment out
				//tmp.push_back(buf.c_str());
				
				if (buf.length() == 0)
					continue;
			    
				if(buf.find("FactorGroup") != string::npos) {  // new factor
					newFactorFound = true;
					idate = 0;
					continue;
				}

				boost::algorithm::split(values, buf, is_any_of(",")); 
				if(newFactorFound) {
					itvec = std::find(h_names.begin(), h_names.end(), fullName);
					if(itvec != h_names.end()) { // set up previous factor
						boost::shared_ptr<Factor> ptr(new Factor); 
						ptr->setAttributes(cur_2_reg_map, "hr_spot", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
						ptr->dataCleanUp(params);
						factors.push_back(std::move(ptr));
						
						// Todo: Comment out
						/*
						ln = tmp.size();
						if(iFactor > 0)
							ln -= 2;
						for(unsigned int k = 0; k < ln; ++k)
							f_clean << tmp[k] << std::endl;
						*/
					}
					++iFactor;
					
					// Todo: Comment out
					//if(tmp.size() > 1) 
					//	tmp.erase(tmp.begin(), tmp.end() - 2);
					
					newFactorFound = false;
					factorGroup = values[0].c_str();
					fullName = values[1].c_str();
					if(fullName == "C00036020_NASDAQ_AAONSpot" || fullName == "C00036110_NYSE_AIRSpot") {
						int tt = 0;
					}
					curr = values[2].c_str();
					hr_dates.clear();
					histDataMap.clear();
					histTermsMap.clear();
				}
				
				if(values.size() == 6) {
					boost::algorithm::split(dayMonYr, values[4], is_any_of("/")); 
					year =atoi(dayMonYr[0].c_str());
					month =atoi(dayMonYr[1].c_str());
					day =atoi(dayMonYr[2].c_str());
					dd = date(year, greg_month(month), day);

					val = atof(values[5].c_str());
						hr_dates.push_back(dd);
						prc.push_back(val);
						histDataMap.insert(std::pair<date, std::vector<double> >(dd, prc));
						histTermsMap.insert(std::pair<date, std::vector<unsigned int> >(dd, term0));
						prc.clear();

					++idate;
				}
			} // end while
			itvec = std::find(h_names.begin(), h_names.end(), fullName);
			if(itvec != h_names.end()) {
				boost::shared_ptr<Factor> ptr(new Factor); 
				ptr->setAttributes(cur_2_reg_map, "hr_spot", factorGroup, fullName, curr, hr_dates, histDataMap, histTermsMap);
				ptr->dataCleanUp(params);
				factors.push_back(std::move(ptr));
				
				// Todo: Comment out
				//for(unsigned int k = 0; k < tmp.size(); ++k)
				//	f_clean << tmp[k] << std::endl;
			}
			// Todo: Comment out
			//tmp.clear();
		}
		inp_file.close();
		// Todo: Comment out
		//f_clean.close();
	}

	return 0;
}
void readShiftInfo(std::string fn, std::map<unsigned int, ShiftScenInfo>& mp) {
	std::string buf(""), name1(""), name2("");;
	std::vector<string> values;
	ifstream ff;
	unsigned int type, scenNum, tenor1, tenor2;
	ShiftType shiftType;
	double shift(0);
	ff.open(fn);
	getline(ff, buf);
	while(getline(ff, buf)) {
		boost::algorithm::split(values, buf, is_any_of(",")); 
		
		if(values.size() < 7) 
			continue;
		
		type = atoi(values[0].c_str());
		scenNum = atoi(values[1].c_str());
		name1 = values[2].c_str();
		name2 = values[3].c_str();
		tenor1 = atoi(values[4].c_str());
		tenor2 = atoi(values[5].c_str());
		shift = atof(values[6].c_str());
		switch(type) {
		case 0: // Base
			shiftType = Base;
			break;
		case 1:
			shiftType = DeltaDn;
			break;
		case 2:
			shiftType = DeltaUp;
			break;
		case 3:
			shiftType = CrossGammaUpUp;
			break;
		case 4:
			shiftType = CrossGammaUpDn;
			break;
		case 5:
			shiftType = CrossGammaDnUp;
			break;
		case 6:
			shiftType = CrossGammaDnDn;
			break;
		}
		ShiftScenInfo inf(name1, name2, tenor1, tenor2, shiftType, shift);
		mp.insert(std::make_pair<unsigned int, ShiftScenInfo> (scenNum, inf));
	}
	
}
bool readVals(std::string fn, std::map<unsigned int, double>& mp) {
	std::ifstream ff;
	std::string buf("");
	std::vector<std::string> values;
	unsigned int scenNum, ii;
	double val;

	ff.open(fn);
	if(!ff.is_open())
		return false;

	getline(ff, buf);
	boost::algorithm::split(values, buf, is_any_of(","));
	if(values.size() < 3)
		return false;

	scenNum = atoi(values[1].c_str());
	while(getline(ff, buf)) {
		boost::algorithm::split(values, buf, is_any_of(","));
		if(values.size() < 3)
			return false;
		ii = atoi(values[1].c_str());
		val = atof(values[2].c_str());
		mp.insert(std::make_pair<unsigned int, double>(ii, val));
	}
	ff.close();
	return true;
}
bool readVals_all(std::string fn, std::map<std::string, std::map<unsigned int, double> >& mp_all, std::vector<std::string>& list) {
	std::ifstream ff;
	std::string buf(""), pfl_name(""), pfl_name_new("");
	std::vector<std::string> values;
	std::map<unsigned int, double> mp;
	unsigned int ii;
	double val;

	ff.open(fn);
	if(!ff.is_open())
		return false;

	getline(ff, buf);

	while(getline(ff, buf)) {
		boost::algorithm::split(values, buf, is_any_of(" ,"));
		if(values.size() < 3)
			return false;

		pfl_name_new = values[0].c_str();
		if(pfl_name_new != pfl_name) {
			if(find(list.begin(), list.end(), pfl_name) != list.end()) {
				mp_all.insert(std::make_pair<std::string, std::map<unsigned int, double> >(pfl_name, mp));
			}
			pfl_name = pfl_name_new;
			mp.clear();
		}
		ii = atoi(values[1].c_str());
		val = atof(values[2].c_str());
		mp.insert(std::make_pair<unsigned int, double>(ii, val));
	}
	ff.close();
	if(find(list.begin(), list.end(), pfl_name) != list.end())
		mp_all.insert(std::make_pair<std::string, std::map<unsigned int, double> >(pfl_name, mp));
	return true;
}
void output_valuations(std::string fn, const std::map<std::string, std::map<unsigned int, double> >& mp_all) {
	std::map<std::string, std::map<unsigned int, double> >::const_iterator itmap;
	std::map<unsigned int, double> mp;
	std::map<unsigned int, double>::iterator it;
	std::string pfl;
	double val;

	std::stringstream out_ff;
	std::ofstream ff;
	ff.open(fn);
	if(!ff.is_open())
		return;

	out_ff.precision(18);
	out_ff << "pfl, scenNum, val" << std::endl;
	for(itmap = mp_all.begin(); itmap != mp_all.end(); ++itmap) {
		pfl = (*itmap).first;
		mp = (*itmap).second;
		for(it = mp.begin(); it != mp.end(); ++it) {
			val = (*it).second;
			out_ff << pfl << "," << (*it).first << "," ;
			//out_ff.precision(18);
			out_ff << val << std::endl;
		}
	}
	
	ff.precision(18);
	ff << out_ff.rdbuf();
	ff.close();
	
	ff.close();
}
void read_pfl_list(std::string fn, std::vector<std::string>& list) {
	std::ifstream ff;
	std::string buf("");
	ff.open(fn);

	if(!ff.is_open())
		return;

	while(getline(ff, buf) ) {
		if(buf != "")
			list.push_back(buf);
	}
	ff.close();
}
// Test CCAR
	
	//ccar_info.read_info(dir_CCAR);
	/*
	FactSet pfl_ccar;
	boost::shared_ptr<Factor> ptr1(new Factor(":C00000117_Index_SPALN_Return")), ptr2(new Factor(":C26099405_Index_INDUA_Return"));
	ptr1->setRegBasedOnCur(cur_2_reg_map);
	ptr1->setType(Eqt_Spot);
	ptr2->setRegBasedOnCur(cur_2_reg_map);
	ptr2->setType(Eqt_Spot);

	//CoreFactor crf(*ptr2);
	std::vector<boost::shared_ptr<Factor> > test_pfl;
	test_pfl.push_back(std::move(ptr1));
	test_pfl.push_back(std::move(ptr2));

	pfl_ccar.addPositions(test_pfl);
	n_proj = ccar_info.getNumProjections();

	//boost::shared_ptr<CoreFactor> crf_levels, crf_vols;
	pfl_ccar.digest_crf_info(ccar_info);

	for(ii = 0; ii < n_proj; ++ii) {
		pfl_ccar.apply_CCAR_info(ccar_info, ii, cur_2_reg_map);
		pfl_ccar.clear_CCAR_info(ccar_info, ii);
	}
	return 1;
	*/
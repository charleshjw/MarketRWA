#include "stdafx.h"

std::map<std::string, std::string> IrCapsFactorsReverseMap = map_list_of 
	(":index_cap1m11", "I07IRVol") (":index_cap1m21","I08IRVol") (":index_cap3m11", "I09IRVol") (":index_cap3m21", "I10IRVol") (":index_cap6m11","I11IRVol") (":index_cap6m21", "I12IRVol") (":index_euro1","I02IRVol");

std::map<std::string, std::string> FxVolFactorsReverseMap=map_list_of 
	(":AUD_USD-Simulation","V03Vol") (":AUD_JPY-Simulation","V02Vol") (":CAD_GBP-Simulation","V22Vol") (":CAD_JPY-Simulation","V04Vol") (":CAD_USD-Simulation","V44Vol") 
	(":CHF_GBP-Simulation","V23Vol") (":CHF_JPY-Simulation","V08Vol") (":CHF_USD-Simulation","V45Vol") (":DKK_GBP-Simulation","V25Vol") (":DKK_USD-Simulation","V47Vol") 
	(":GBP_JPY-Simulation","V29Vol") (":GBP_NOK-Simulation","V31Vol") (":GBP_USD-Simulation","V33Vol") (":HKD_USD-Simulation","V52Vol") (":JPY_USD-Simulation","V55Vol") 
	(":MXN_USD-Simulation","V56Vol") (":NOK_USD-Simulation","V58Vol") (":NZD_USD-Simulation","V39Vol") (":RUB_USD-Simulation","V60Vol") (":SEK_USD-Simulation","V62Vol") 
	(":USD_ZAR-Simulation","V63Vol") (":CHF_SEK-Simulation","V09Vol") (":SEK_DKK-Simulation","V40Vol") (":EUR_AUD-Simulation","V73Vol") (":EUR_CHF-Simulation","V71Vol") 
	(":EUR_GBP-Simulation","V67Vol") (":EUR_JPY-Simulation","V68Vol") (":EUR_NOK-Simulation","V69Vol") (":EUR_SEK-Simulation","V70Vol") (":EUR_USD-Simulation","V66Vol") 
	(":EUR_CZK-Simulation","V72Vol") (":EUR_TRY-Simulation","V74Vol") (":USD_TRY-Simulation","V75Vol");



void test_hist_vs_simulated(std::string dir, std::vector<boost::shared_ptr<Factor> >& facs) {
	unsigned int i, j, nsim(facs[0]->_sim_rets.size1()), nrets(facs[0]->_rets.size1()), n1, n2;
	std::vector<boost::shared_ptr<Factor> >::iterator itfac, jtfac;
	std::ofstream ff;
	std::string fn(dir + "test_correl.csv");
    ff.open(fn.c_str());

	matrix<double> h_rets1, h_rets2, s_rets1, s_rets2, h_cov, h_cor, s_cov, s_cor, tmp, h_cov1, h_cov2, s_cov1, s_cov2;
	boost::numeric::ublas::vector<double> h_vol1_inv, h_vol2_inv, s_vol1_inv, s_vol2_inv;

	ff << "Factor1,Factor2,Name1, Value1, Name2, Value2, Name3, Value3" << std::endl;
	for(n1 = 0; n1 < facs.size(); ++n1) {
		h_rets1 = facs[n1]->_rets;
		s_rets1 = facs[n1]->_sim_rets;
		h_cov1 = prod(trans(h_rets1), h_rets1) / nrets;
		s_cov1 = prod(trans(s_rets1), s_rets1) / nsim;

		h_vol1_inv.clear();
		h_vol1_inv.resize(h_cov1.size1());
		for(i = 0; i < h_cov1.size1(); ++i)
			h_vol1_inv[i] = 1.0/sqrt(h_cov1(i,i));
		s_vol1_inv.clear();
		s_vol1_inv.resize(s_cov1.size1());
		for(i = 0; i < s_cov1.size1(); ++i)
			s_vol1_inv[i] = 1.0/sqrt(s_cov1(i,i));

		for(n2 = n1; n2 < facs.size(); ++n2) {
			h_rets2 = facs[n2]->_rets;
			s_rets2 = facs[n2]->_sim_rets;
			h_cov2 = prod(trans(h_rets2), h_rets2) / nrets;
			s_cov2 = prod(trans(s_rets2), s_rets2) / nsim;

			h_vol2_inv.clear();
			h_vol2_inv.resize(h_cov2.size1());
			for(i = 0; i < h_cov2.size1(); ++i)
				h_vol2_inv[i] = 1.0/sqrt(h_cov2(i,i));
			s_vol2_inv.clear();
			s_vol2_inv.resize(s_cov2.size1());
			for(i = 0; i < s_cov2.size1(); ++i)
				s_vol2_inv[i] = 1.0/sqrt(s_cov2(i,i));

			h_cov = prod(trans(h_rets1), h_rets2) / nrets;
			s_cov = prod(trans(s_rets1), s_rets2) / nsim;

			//banded_matrix<double> s_var(s_cov.size1(), s_cov.size1(), 0, 0), h_var(h_cov.size1(), h_cov.size1(), 0, 0);
			h_cor.clear();
			h_cor.resize(h_cov.size1(), h_cov.size2());
			for(i = 0; i < h_cov.size1(); ++i)
				for(j = 0; j < h_cov.size2(); ++j)
					h_cor(i, j) = h_cov(i, j) * h_vol1_inv[i] * h_vol2_inv[j];

			s_cor.clear();
			s_cor.resize(s_cov.size1(), s_cov.size2());
			for(i = 0; i < s_cov.size1(); ++i)
				for(j = 0; j < s_cov.size2(); ++j)
					s_cor(i, j) = s_cov(i, j) * s_vol1_inv[i] * s_vol2_inv[j];
			
			ff << facs[n1]->getHistName() << "," << facs[n2]->getHistName() << ",Hist Vol1,";
			for(i = 0; i < h_vol1_inv.size(); ++i)
				ff << 1.0 / h_vol1_inv[i] << ",";
			ff << "Hist Vol2,";
			for(i = 0; i < h_vol2_inv.size(); ++i)
				ff << 1.0 / h_vol2_inv[i] << ",";
			ff << "Hist Cor,";
			print_matrix(ff, h_cor);
			ff << std::endl;
			ff << facs[n1]->getHistName() << "," << facs[n2]->getHistName() << ",Sim Vol1,";
			for(i = 0; i < s_vol1_inv.size(); ++i)
				ff << 1.0 / s_vol1_inv[i] << ",";
			ff << "Sim Vol2,";
			for(i = 0; i < s_vol2_inv.size(); ++i)
				ff << 1.0 / s_vol2_inv[i] << ",";
			ff << "Sim Cor,";
			print_matrix(ff, s_cor);
			ff << std::endl;

			ff << facs[n1]->getHistName() << "," << facs[n2]->getHistName() << ",Dif Vol1,";
			for(i = 0; i < s_vol1_inv.size(); ++i)
				ff << s_vol1_inv[i] / h_vol1_inv[i] -1 << ",";
			ff << "Dif Vol 2,";
			for(i = 0; i < s_vol2_inv.size(); ++i)
				ff << s_vol2_inv[i] / h_vol2_inv[i] -1 << ",";
			ff << "Dif Corr,";

			tmp = s_cor - h_cor;
			print_matrix(ff, tmp);
			ff << std::endl;
		}
	}
	ff.close();
}

void test_hist_vs_simulated_Mapped(std::string dir, std::vector<boost::shared_ptr<Factor> >& facs1, std::vector<boost::shared_ptr<Factor> >& facs2, unsigned int sim_horizon) {
	unsigned int i, j, nsim(facs1[0]->_sim_rets.size1()), nrets, n1, n2;
	std::vector<boost::shared_ptr<Factor> >::iterator itfac, jtfac;
	std::ofstream ff;
	std::string fn(dir + "test_correl_mapped.csv");
    ff.open(fn.c_str());

	matrix<double> h_rets1, h_rets2, s_rets1, s_rets2, h_cov, h_cor, s_cov, s_cor, tmp, h_cov1, h_cov2, s_cov1, s_cov2;
	boost::numeric::ublas::vector<double> h_vol1_inv, h_vol2_inv, s_vol1_inv, s_vol2_inv;
	std::vector<date> dts1, dts2, dts_res, new_dts;

	ff << "Factor1,Factor2,Name1, Value1, Name2, Value2, Name3, Value3" << std::endl;
	for(n1 = 0; n1 < facs1.size(); ++n1) {
		s_rets1 = facs1[n1]->_sim_rets;
		s_cov1 = prod(trans(s_rets1), s_rets1) / nsim;
		s_vol1_inv.clear();
		s_vol1_inv.resize(s_cov1.size1());
		for(i = 0; i < s_cov1.size1(); ++i)
			s_vol1_inv[i] = 1.0/sqrt(s_cov1(i,i));

		dts1 = facs1[n1]->getHrDates();

		for(n2 = n1; n2 < facs2.size(); ++n2) {
			s_rets2 = facs2[n2]->_sim_rets;
			s_cov2 = prod(trans(s_rets2), s_rets2) / nsim;s_cov2 = prod(trans(s_rets2), s_rets2) / nsim;
			s_vol2_inv.clear();
			s_vol2_inv.resize(s_cov2.size1());
			for(i = 0; i < s_cov2.size1(); ++i)
				s_vol2_inv[i] = 1.0/sqrt(s_cov2(i,i));

			s_cov = prod(trans(s_rets1), s_rets2) / nsim;
			s_cor.clear();
			s_cor.resize(s_cov.size1(), s_cov.size2());
			for(i = 0; i < s_cov.size1(); ++i)
				for(j = 0; j < s_cov.size2(); ++j)
					s_cor(i, j) = s_cov(i, j) * s_vol1_inv[i] * s_vol2_inv[j];

			// Find common dates and returns for these dates
			dts2 = facs2[n2]->getHrDates();
			std::vector<date>::iterator it_res;
			
			dts_res.resize(max(dts1.size(), dts2.size()));
			it_res = set_intersection(dts1.begin(), dts1.end(), dts2.begin(), dts2.end(), dts_res.begin());
			dts_res.resize(it_res - dts_res.begin());
			
			facs1[n1]->setRetsOnDates(dts_res, sim_horizon, h_rets1, new_dts);
			facs2[n2]->setRetsOnDates(dts_res, sim_horizon, h_rets2, new_dts);
			nrets = dts_res.size() - 1;

			h_cov1 = prod(trans(h_rets1), h_rets1) / nrets;
			h_vol1_inv.clear();
			h_vol1_inv.resize(h_cov1.size1());
			for(i = 0; i < h_cov1.size1(); ++i)
				h_vol1_inv[i] = 1.0/sqrt(h_cov1(i,i));
				
			
			h_rets2 = facs2[n2]->_rets;
			h_cov2 = prod(trans(h_rets2), h_rets2) / nrets;
			h_vol2_inv.clear();
			h_vol2_inv.resize(h_cov2.size1());
			for(i = 0; i < h_cov2.size1(); ++i)
				h_vol2_inv[i] = 1.0/sqrt(h_cov2(i,i));

			h_cov = prod(trans(h_rets1), h_rets2) / nrets;

			//banded_matrix<double> s_var(s_cov.size1(), s_cov.size1(), 0, 0), h_var(h_cov.size1(), h_cov.size1(), 0, 0);
			h_cor.clear();
			h_cor.resize(h_cov.size1(), h_cov.size2());
			for(i = 0; i < h_cov.size1(); ++i)
				for(j = 0; j < h_cov.size2(); ++j)
					h_cor(i, j) = h_cov(i, j) * h_vol1_inv[i] * h_vol2_inv[j];

			ff << facs1[n1]->getHistName() << "," << facs2[n2]->getHistName() << "," << nrets << ",Hist Vol1,";
			for(i = 0; i < h_vol1_inv.size(); ++i)
				ff << 1.0 / h_vol1_inv[i] << ",";
			ff << "Hist Vol2,";
			for(i = 0; i < h_vol2_inv.size(); ++i)
				ff << 1.0 / h_vol2_inv[i] << ",";
			ff << "Hist Cor,";
			print_matrix(ff, h_cor);
			ff << std::endl;
			ff << facs1[n1]->getHistName() << "," << facs2[n2]->getHistName() << "," << nrets << ",Sim Vol1,";
			for(i = 0; i < s_vol1_inv.size(); ++i)
				ff << 1.0 / s_vol1_inv[i] << ",";
			ff << "Sim Vol2,";
			for(i = 0; i < s_vol2_inv.size(); ++i)
				ff << 1.0 / s_vol2_inv[i] << ",";
			ff << "Sim Cor,";
			print_matrix(ff, s_cor);
			ff << std::endl;

			ff << facs1[n1]->getHistName() << "," << facs2[n2]->getHistName() << ",Dif Vol1,";
			for(i = 0; i < s_vol1_inv.size(); ++i)
				ff << s_vol1_inv[i] / h_vol1_inv[i] -1 << ",";
			ff << "Dif Vol 2,";
			for(i = 0; i < s_vol2_inv.size(); ++i)
				ff << s_vol2_inv[i] / h_vol2_inv[i] -1 << ",";
			ff << "Dif Corr,";

			tmp = s_cor - h_cor;
			print_matrix(ff, tmp);
			ff << std::endl;
		}
	}
	ff.close();
}
void print_matrix(std::ofstream& ff, matrix<double>& mat) {
	unsigned int n1(mat.size1()), n2(mat.size2()), i1, i2;

	for(i1 = 0; i1 < n1; ++i1)
		for(i2 = 0; i2 < n2; ++i2)
			ff << mat(i1, i2) << ",";
}

int output_scenarios(std::string dir, std::string asOfDate, unsigned int nsim, std::vector<boost::shared_ptr<Factor> >& facs, bool output_diff) {

	std::vector<unsigned int> terms;
	std::vector<double> eod_val;
	matrix<double> rets;
	std::ofstream out_file;
	double prob(1.0/nsim);
	std::string s_name;
	unsigned int i, j, k, nterms;
	std::ostringstream ostr;
		
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	std::string scen_color("green"), scenSetName("msmc 2000"), scen_name("Mmc2000"), TimeEvolution("Constant"), TimeEvolutionToTrigger("Constant"),
		TriggerHolder("@Standard()"), ScenShiftRule("Trigger Time"), ScenType("non-parallel shift"), ScenReplVal("Term");

	std::vector<boost::shared_ptr<Factor> >::iterator it, jt;
	std::string fn(dir + "scen_MC.csv");
    out_file.open(fn.c_str());
	
	out_file << "Scenario Set,Scenario Set Name,Scenario Name,Scenario Probability,Scenario Color,Scenario Variable,Scenario Start Time,Scenario Attribute,Time Evolution,Time Evolution to Trigger,Trigger Holder,Scenario Shift Rule,Scenario Type,Scenario Replacement Value" << std::endl;

	for(i = 0; i < nsim; ++i) {
		j = 0;
		std::ostringstream ostr;
		ostr << (i + 1);
		for(itfac = facs.begin(); itfac != facs.end(); ++itfac) {
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
						if((*itfac)->isLogNormal())
							out_file << ",,,,,,,,,,,,," << terms[k] << "," << eod_val[k] * (exp(rets(i,k)) - 1.0) << std::endl;
						else
							out_file << ",,,,,,,,,,,,," << terms[k] << "," << rets(i,k) << std::endl;
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
	}
	// Base scenario
	j = 0;
	std::ostringstream ostr2;
	ostr2 << (i + 1);
	for(itfac = facs.begin(); itfac != facs.end(); ++itfac) {
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

 	return 0;
}
std::vector<boost::shared_ptr<Factor> >::iterator findFactor(FactType iType, std::string subType, std::string cur, std::vector<boost::shared_ptr<Factor> >& facts) {

	std::vector<boost::shared_ptr<Factor> >::iterator itfac(facts.begin());
	while(itfac != facts.end() && !((*itfac)->getType() == iType && (*itfac)->getCurr() == cur && (*itfac)->getSubType() == subType))
		++itfac;

	return itfac;
}


void readStockmapFile(std::string stockmap_fn, std::map<unsigned int, std::string>& stockMap) {
	std::ifstream stockmap_file;
	stockmap_file.open(stockmap_fn);
	unsigned int ind(0);
	std::string buf(""), h_name(""), saved_name("");
	std::vector<std::string> values, sub_values;
	using boost::is_any_of;

	while(getline(stockmap_file, buf)) {
		// File structure: "ACS_C00819010	0	USD	1" - equityName_cusip interger currency tenor
		boost::algorithm::split(values, buf, is_any_of("|"));	
		
		if(values.size() == 4 && values[0] != saved_name) {
			boost::algorithm::split(sub_values, values[0], is_any_of("_"));
			ind = atoi(values[1].c_str());
			if(sub_values.size() > 2) {
				h_name = values[0] + "Spot";
				stockMap.insert(std::pair<unsigned int, std::string>(ind, h_name));
				saved_name = values[0];
			}
		}
	}
}


std::string s_name_2_h_name(std::string s_name) {
	 // Convert an "s_name" to "h_name"
	 // This should be done in robust way!!! Tanya
	std::vector<std::string> values, split;
	std::string s1, res, subType;
	unsigned int n;
	std::map<std::string, std::string>::iterator itmap;

	res = "ZZZ_" + s_name; // Good way to find not handled names Tanya
	//res = "";
	boost::algorithm::split(values, s_name, is_any_of(":"));

	if(s_name == ":FXEUR")
		int tt = 0;

	if(values.size() == 2) {
		s1 = values[1];
		if(s1 == "REFI_DAILY" || s1 =="TURNOVER_DAILY" || s1 == "CDS_RADIOSHA_7C547B_SNRFOR_MR_EQTV") {
			res = s1;
			return res;
		}
		boost::algorithm::split(values, s1, is_any_of("_"));
		n = values.size();
		
		switch (n) {
		case 1:																// FX Spot					//	:FXEUR						EURSpot
			if(values[0].substr(0,2)=="FX")					
				res = values[0].substr(2) + "Spot";
			break;
		case 2:	
			if(values[0] =="factorGB") { // :factorGB_0		J00IRVol
				int nn = atoi(values[1].c_str());
				std::ostringstream ostr;
				ostr << nn;
				if(nn < 10)
					res = std::string("J0" + ostr.str() + "IRVol");
				else
					res = std::string("J" + ostr.str() + "IRVol");
			}
			else if(values[1] == "Forward" && values[0].substr(0,2) == "FX")		// FX Forward				// :FXEUR_Forward				EURForward
				res = values[0].substr(2) +  values[1];
			else {	
				if(values[0] == "index") {									// Cap volatility			// :index_cap1m11				I07IRVol
					itmap = IrCapsFactorsReverseMap.find(s_name);
					if(itmap != IrCapsFactorsReverseMap.end())
						res = itmap->second;
				}
				else {
					boost::algorithm::split(split, values[1], is_any_of("-"));
					if(split.size()== 2 && split[1] == "Simulation")  {		// FX pair volatility		// :JPY_USD-Simulation			V55Vol
						itmap = FxVolFactorsReverseMap.find(s_name);
						if(itmap != FxVolFactorsReverseMap.end())
							res=itmap->second;
						else {// try to find if :GBP_USD-Simulation exists instead of :USD_GBP-Simulation
							std::string cur1(split[0]), cur2(values[0]);
							std::string s_name2(":" + cur1 + "_" + cur2 + "-Simulation");
							itmap = FxVolFactorsReverseMap.find(s_name2);
							if(itmap != FxVolFactorsReverseMap.end())
								res=itmap->second;
						}
					}
					else if(values[0].substr(0,2) == "IR")	{				// IR curves				// :IREUR_Interbank				EURIntebank	
						subType = values[1];
						if(subType == "Mmarket")
							subType = "MM";
						res = values[0].substr(2) + subType;
					}
				}
			}
			break;
		case 3:																							
			if(values[2] == "EQTV")										// Equity volatilities			// :AAPL_C03783310_EQTV			AAPL_C03783310_EQTV
				res = s1;
			break;
		case 4:
			if(s1 == "IRUSD_MBS_CC_OAS")
				res = "USDMBS_CC_OAS";
			else if(values[3] == "Return") { // h_name == "C85299600_OTC_SPTRSpot", _s_name = ":C85299600_Index_SPTR_Return"; This should be redone
				if(s1 == "C85299600_Index_SPTR_Return")
					res = "C85299600_OTC_SPTRSpot";
				else if(s1 == "C12496G10_Index_RUSSELL3000_Return")
					res = "C12496G10_OTC_RUSSELL3000Spot";
				else if (s1 == "C00000117_Index_SP500_Return")
					res = "C00000117_Index_SPALNSpot";
				else if(s1 == "C12497K10_Index_VIX_Return")
					res = "C12497K10_OTC_VIXSpot";
				else if(s1 == "C26099405_Index_INDUA_Return")
					res = "C26099405_OTC_INDUASpot";
				else
					res = values[0] + "_" + values[1] + "_" + values[2] + "Spot";
			}
			break;
		case 5:
		case 6:
		case 7:
		case 8:
			if(s1.length() > 7 && s1.substr(s1.length() - 7, 7) == "_Return") 
				res = s1.substr(0, s1.length()- 7) + "Spot";
			//if(values[4] == "Return") // _h_name = "C78462F10_NYSE_SPY_ETFSpot", _s_name = ":C78462F10_NYSE_SPY_ETF_Return"
			//	res = values[0] + "_" + values[1] + "_" + values[2] + "_" + values[3] + "Spot";
			break;
		}
	}

	return res;
}

void readEmpirDistFactList(std::string fn, std::map<std::string, std::vector<unsigned int> >& empirDistMap) {
	std::ifstream ff;
	ff.open(fn);
	std::string buf, name;
	std::vector<std::string> values;
	unsigned int tenor;
	std::map<std::string, std::vector<unsigned int> >::iterator itmap;
	std::vector<unsigned int> tenors;

	while(getline(ff, buf)) {
		if(buf.find("#") != std::string::npos) 
			continue;

		boost::algorithm::split(values, buf, is_any_of(","));
		if(values.size() == 2) {
			name = values[0];
			tenor = atoi(values[1].c_str());
			itmap = empirDistMap.find(name);
			std::vector<unsigned int> updated_tenors;
			if(itmap == empirDistMap.end()) {
				updated_tenors.push_back(tenor);
				//empirDistMap.insert(std::make_pair<std::string, unsigned int>(name, tenor));
			}
			else {
				tenors = itmap->second;
				updated_tenors = tenors;
				updated_tenors.push_back(tenor);
			}
			empirDistMap[name] = updated_tenors;	
		}
	}
	ff.close();
}
std::vector<boost::shared_ptr<Factor> >::const_iterator matchNameWithFactor(const std::string name, const std::vector<boost::shared_ptr<Factor> >& facts, bool isHistName) {

	std::vector<boost::shared_ptr<Factor> >::const_iterator itfac(facts.begin());
	if(isHistName)
		while(itfac != facts.end() && name != (*itfac)->getHistName())
			++itfac;
	else
		while(itfac != facts.end() && name != (*itfac)->getScenName())
			++itfac;
	return itfac;
}
/*
void readEOD(std::vector<boost::shared_ptr<Factor> >& facs, std::string dir) {

	std::ifstream ff;
	std::string fn[]={"ir_zerocurves.csv", "fx_zerocurves.csv", "other_zerocurves.csv", "creditbond_zerocurves.csv", "IndexCurves.csv"};  //Temp Tanya
	static const std::vector<std::string> file_names(fn, fn + sizeof(fn)/sizeof(fn[0]));
	std::vector<unsigned int> ncolumns(20), terms;
	unsigned int i, iterm, numFiles(file_names.size());
	std::string buf(""), s_name("");
	std::vector<std::string> values;
	bool newFactorFound(false), factorFound(false);
	std::vector<boost::shared_ptr<Factor> >::iterator itfac(facs.end()), jtfac;
	std::vector<double> arg, vals;
	std::map<double, double> mp;
	std::map<double, double>::iterator it, jt;
	double vv(0.0);

	// Temp Tanya this should be changes
	ncolumns[0] = 13; // is_zerocurve.csv
	ncolumns[1] = 12; // fx_zerocurves.csv
	ncolumns[2] = 12; // other_zerocurves.csv
	ncolumns[3] = 12; // creditbond_zerocurves.csv
	ncolumns[4] = 12; // IndexCurves.csv

	for(i = 0; i < numFiles; ++i) {
		ff.open(dir + file_names[i]);
		
		while(getline(ff, buf)) {
			if(buf.find("RiskMetrics Link") != std::string::npos) {
				if(itfac != facs.end()) { // Process the factor
					if((*itfac)->getHistName() == "C78462F10_NYSE_SPY_ETFSpot")
						int tt = 0;
				
					if((*itfac)->getUseProxyEOD()) {
						jtfac = (*itfac)->findProxyFact(facs);
						if((*itfac)->getType() == (*jtfac)->getType())
							(*itfac)->setEODVal((*jtfac)->getEODVal());
					}
					else {
						terms = (*itfac)->getTerms();
						arg.resize(terms.size());
						vals.resize(terms.size());
						for(iterm = 0; iterm < terms.size(); ++iterm)
							arg[iterm] = terms[iterm];
						interpolate_fwd_linear(mp, arg, vals);
					
						(*itfac)->setEODVal(vals);
						vals.clear();
						mp.clear();
					}
				}
				newFactorFound = true;
				s_name = "";
				itfac = facs.end();
				continue;
			}
			boost::algorithm::split(values, buf, is_any_of(",")); 
			if(newFactorFound) {
				if(values.size() > 3) {
					s_name = ":" + values[2];
					itfac = matchNameWithFactor(s_name, facs, false);
					newFactorFound = false;
				}
				continue;
			}

			if(itfac != facs.end()) {
				if(values.size() == ncolumns[i]) 
					mp.insert(std::make_pair(atoi(values[ncolumns[i] - 2].c_str()), atof(values[ncolumns[i] - 1].c_str())));
			}
		}
		ff.close();
	}
}
*/
/*
void setEOD_fromHist(std::vector<boost::shared_ptr<Factor> >& facs) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	unsigned int iterm, nterms;
	std::
	for(itfac = facs.begin(); itfac != facs.end(); ++itfac) {
		nterms = (*itfac)->getNterms();

	}
}
(/

/*
void check_factors(std::vector<boost::shared_ptr<Factor> >& facs, std::string fn) {
	// Adhoc measure until reading eod data is implemented this is for testing ONLY!!! Tanya
	double rate_eod(0.01), vol_eod(0.3);
	std::vector<boost::shared_ptr<Factor> >::iterator itfac, jtfac;
	FactType type;
	unsigned int i, nterms;
	std::vector<double> eod, vol;
	std::vector<unsigned int> terms;
	std::string proxy_name("");
	std::ofstream ff;

	ff.open(fn);
	ff << "h_name,  Proxy name, mapping_type, how eod was assigned, eod val, eod size " << std::endl;
	for(itfac = facs.begin(); itfac != facs.end(); ++itfac) {
		
		if((*itfac)->getHistName() == "USDLIB1Y")
			int tt = 0;

		nterms = (*itfac)->getNterms();
		terms = (*itfac)->getTerms();
		
		else if((*itfac)->getEODVal().size() != nterms){
			ff << (*itfac)->getHistName() << "," << (*itfac)->getProxyName() << "," << (*itfac)->getType() << ",EOD size is differenr from nterms," << nterms << "," << eod.size() << std::endl;
		}
	}
	ff.close();
}
*/
/* Test
	boost::mt19937                     gener(1);
    boost::normal_distribution<double> normal_dist;   // Normal Distribution
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal_dist);
	rng.engine().seed(1);
	unsigned int n1(500), n2(1000);
	std::vector<double> empir(n1), vnorm(n2), res(n2);
	double mean, std;
	for(i = 0; i < n1; ++i) 
		empir[i] = 100 * rng() + 5;
	for(i = 0; i < n2; ++i) 
		vnorm[i] = 90 * rng() + 4;

	transform_2_empirical(empir, vnorm, res);
	normalize_vec(res, mean, std);
	normalize_vec(res, mean, std);
	*/
/*
	std::string check_fn(dir_output + "check_eod_val.csv");
	check_factors(h_factors, check_fn);
	*/
	// Print factors
	/*
	cleaned_facs.open(cleaned_facs_fn);
	for(itfac = h_factors.begin(); itfac != h_factors.end(); ++itfac)
		(*itfac)->print_hist(cleaned_facs);
	cleaned_facs.close();
	*/
/*
void getFactorsByTypeAndDataAval(std::vector<boost::shared_ptr<Factor> >& factors, HistDataAvailabilityType dataAvalType, FactType factType, std::vector<boost::shared_ptr<Factor> >& f_vec) {
	std::vector<boost::shared_ptr<Factor> >::iterator itfac;
	for(itfac = factors.begin(); itfac != factors.end(); ++itfac) {
		if((*itfac)->getType() == factType && (*itfac)->getDataAvalType() == dataAvalType && (*itfac)->getNterms() > 0)
			f_vec.push_back(*itfac);
	}
			
}
*/
/*
void output_factors_info(std::string dir_output, std::vector<boost::shared_ptr<Factor> >& factors) {
	
	std::ofstream log_file;
	std::string hist_fact_info_fn("hist_fact_info.csv");
	std::string fn(dir_output + hist_fact_info_fn);
    log_file.open(fn.c_str());

	unsigned int i;

	log_file << "Number, Factor name, Mapping type, Factor type, Currency, Region, Number of dates, Number of terms, vol, Resid vol, R2, Regr coef" << std::endl;
	for(i = 0; i < factors.size(); ++i) {
		log_file << i << ",";
		factors[i]->print_factor_info(log_file);
	}
	
	log_file.close();
}
*/
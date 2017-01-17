#ifndef CURRENCYMAPS_H
#define CURRENCYMAPS_H

class CurrencyMaps {
public:
	bool read_region_2_main_curr();
	bool read_curr_to_reg();
	bool validate_based_on_data_availablility();
	bool create_currency_2_factor_map();
	bool getMainFactor(std::string currency);

private:
	std::map<std::string, std::string> region_2_main_curr;
	std::map<std::string, std::string> curr_2_reg;
	std::string region_2_main_curr_fn;	// file name to get map
	std::string curr_2_reg_fn;			// file name to get map
	std::vector<boost::shared_ptr<Factor> > _fxSpotMainFactors;
	std::vector<boost::shared_ptr<Factor> > _fxForwardMainFactors;
	std::vector<boost::shared_ptr<Factor> > _irCurveMainFactors;
};

#endif CURRENCYMAPS_H
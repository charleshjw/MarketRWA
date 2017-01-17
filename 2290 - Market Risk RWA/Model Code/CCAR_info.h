#ifndef CCARINFO_H
#define CCARINFO_H

class CCAR_shift;

class CCAR_info {

public:
	CCAR_info()  : _num_projections(1) {};
	~CCAR_info() {};
	bool			read_info(std::string dir);
	bool			read_reg_2_core_cur_map(std::string fn);
	bool			read_rfTypeCur_coreRf_map(std::string fn);
	bool			read_core_rf_info(std::string fn);
	bool			read_core_rf_projections(std::string fn);

	unsigned int	getNumProjections() const {return _num_projections;};
	//bool			read_proj(std::string fn, std::map<std::string, matrix<double> >& map);
	bool			find_projection(std::string s_name,  bool isEOD, unsigned int i_proj, std::map<double, double>& proj) const;
	//bool			find_mapping(std::string name, bool isEOD, std::string& mapped_name, double& shift, bool& isMult);

	const std::map<std::pair<std::string, std::string>, std::string>&			getFactTypeCur_coreFact_mp() const {return _rfTypeCur_coreRf_mp;};
	const std::map<std::string, std::string>&									getReg_2_cur_mp() const {return _reg_2_cur_mp;};
	const std::map<std::string, bool>&											get_core_rf_info() const {return _core_fact_info;};
	const std::map<std::string, std::map<unsigned int, std::vector<double> > >&	get_eod_proj() const {return _eod_proj;};
	void calc_crf_shifts(const FactSet& facs, std::vector<boost::shared_ptr<CCAR_shift> >& ccar_shifts);
private:
	unsigned int _num_projections;
	std::map<std::string, std::map<unsigned int, std::vector<double> > >  _eod_proj; // for each core factor the pair of tenor and projection vector
	std::map<std::string, std::map<unsigned int, std::vector<double> > >  _vol_proj; // for each core factor the pair of tenor and projection vector
	std::map<std::string, bool> _core_fact_info;
	std::map<std::pair<std::string, std::string>, std::string> _rfTypeCur_coreRf_mp;
	std::map<std::string, std::string> _reg_2_cur_mp;
};

#endif
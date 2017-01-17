#ifndef UTILITIES_H
#define UTILITIES_H

void	test_hist_vs_simulated(std::string dir, std::vector<boost::shared_ptr<Factor> >& factors);
void	test_hist_vs_simulated_Mapped(std::string dir, std::vector<boost::shared_ptr<Factor> >& facs1, std::vector<boost::shared_ptr<Factor> >& facs2, unsigned int sim_horizon);
void	print_matrix(std::ofstream& ff, matrix<double>& mat);
void	readStockmapFile(std::string stockmap_fn, std::map<unsigned int, std::string>& stockMap);
void	getFactorsByType (std::vector<boost::shared_ptr<Factor> >& factors, FactType factType, std::vector<boost::shared_ptr<Factor> >& f_vec);
void	readEmpirDistFactList(std::string fn, std::map<std::string, std::vector<unsigned int> >& empirDistMap);
void	assignEmpirDistFlag(std::vector<boost::shared_ptr<Factor> >& facs, std::map<std::string, std::vector<unsigned int> >& empirDistMap);
//bool	interpolate_linear(std::vector<double>& x, std::vector<double>&y, std::vector<double>& arg, std::vector<double>& res);
std::vector<boost::shared_ptr<Factor> >::const_iterator matchNameWithFactor(const std::string h_name, const std::vector<boost::shared_ptr<Factor> >& facts, bool isHistName);
std::string  s_name_2_h_name(std::string s_name);

#endif

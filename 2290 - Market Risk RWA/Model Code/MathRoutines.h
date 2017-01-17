#ifndef MATHROUTINES_H
#define MATHROUTINES_H

#include <boost/numeric/ublas/lu.hpp> 
#include "boost/random.hpp"
#include <boost/numeric/ublas/vector_proxy.hpp> // To compile lu_factorize
#include <boost/numeric/ublas/triangular.hpp> // To compile lu_factorize

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) > (b) ? (b) : (a))

using boost::math::normal; // typedef provides default type is double.
bool	regress(matrix<double>& X, matrix<double>& Y, matrix<double>& coef, std::vector<double>& vols_resid, std::vector<double>& R2);
void	svdcmp(matrix<double>& a, std::vector<double>& w, matrix<double>& v);
int		genIndepNorm(int m, std::vector<double>& vec);
void	sortWithIndex(std::vector<double>& vec, std::vector<unsigned int>& index) ;
bool	pairsValComp(const std::pair<int, double>& p1, const std::pair<int, double> p2);
double	orthog_check(matrix<double>& randoms);
void	gramm_shmidt_orth(matrix<double>& mat);
double	stdev(std::vector<double>& vec, std::vector<double>& wts);
void	calculate_weights(double decay, std::vector<date>& dts, std::vector<double>& wts, unsigned int sim_horizon);
bool	InvertMatrix (const matrix<double>& input, matrix<double>& inverse);
void	svdLanczos(matrix<double>& histRets, unsigned int nc, matrix<double>& eigenvec, std::vector<double>& eigenval);
double	fit_normal_cdf_std(double x, double y);
double	cdf(double x);
void	transform_2_empirical(std::vector<double>& empir_vec, std::vector<double>& norm_vec, std::vector<double>& res);
void	normalize_vec(std::vector<double>& vec, double& mean, double& vol);
double	interpolate_linear(std::vector<double>& x, std::vector<double>&y, double val);
void	interpolate_linear(std::map<double,double>& mp, std::vector<double>& arg, std::vector<double>& res);
void	interpolate_fwd_linear(std::map<double,double>& mp, std::vector<double>& arg, std::vector<double>& res);
void	interpolate_quadratic(std::vector<double>& xi, std::vector<double>& yi, std::vector<double>& xo, std::vector<double>& yo);
#endif MATHROUTINES
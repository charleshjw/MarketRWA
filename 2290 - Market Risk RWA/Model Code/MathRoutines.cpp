#include "stdafx.h"

static double at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void svdLanczos(matrix<double>& histRets, unsigned int nPC, matrix<double>& eigenvec_out, std::vector<double>& eigenval_out)
{
	using namespace boost::numeric::ublas;
	// Historical returns matrix histRets, np (number of factors) rows, nd (number of returns) columns
	unsigned int i, j, icount(0), k, nc(nPC * 2), nEig(nc), nSelect(nc), np(histRets.size2()), nd(histRets.size1());
	double norm, eps(1E-11), acc;
	boost::numeric::ublas::vector<double> vols(np), q0(np), q1(np), q2(np), vec(nd), al(nc), bt(nc+1), wts(nc + 1);
	std::vector<double> rnd, eigenval;
	matrix<double> QQ(np, nc), tridiag(nc, nc), ww(nc+1, nc+1), rets(nd, np), eigenvec;

	eigenval.resize(nc);
	eigenvec.resize(nc, nc);

	// Scale historical returns with volatility
	for(i = 0; i < np; ++i) {
		vols[i] = sqrt(inner_product(column(histRets,i).begin(), column(histRets,i).end(), column(histRets,i).begin(), 0.0));
		for(k =0; k < nd; ++k)
			rets(k, i) = histRets(k, i) / vols[i];
	}

	// Set q1
	genIndepNorm(np, rnd);
	for(i = 0; i < np; ++i)
		q1[i] = rnd[i];

	bt[0] = 0.0;
	norm = inner_product(q1.begin(), q1.end(), q1.begin(), 0.0);
	norm = 1.0/sqrt(norm);
	q1 *= norm;
	ww(0, 0) = 1.0;
	
	// Lanczos tri-diagonalization
	acc = 0.0;
	for(i=0; i < nc; ++i) {
		column(QQ, i) = q1;

		// Lanczos step: Aq(i) = bt(i) * q(i-1) + a(i) * q(i) + bt(i+1) * q(i+1)
		vec = prod(rets, q1);
		q2 = prod(trans(rets), vec);
		
		al[i] = inner_product(q1.begin(), q1.end(),q2.begin(), 0.0);
		acc += al[i];
		q2 -= al[i] * q1;
		if(i > 1)
			q2 -= bt[i] * q0;

		bt[i+1] = inner_product(q2.begin(), q2.end(), q2.begin(), 0.0);
		bt[i+1] = sqrt(bt[i+1]);
		if(abs(bt[i+1]) < eps)
			break;

		q2 /= bt[i+1];
		
		// Partial reorthogonalization: weights computation
		ww(i+1,i+1) = 1.0;
		
		// Update weights
		wts = prod(trans(QQ), q2);
		for(k = 0; k <= i; ++k) {
			ww(i+1,k)= wts[k];
			ww(k,i+1) = wts[k];
		}
		// Update vectors
		for(k = 0; k <= i; ++k) 
			q2 -= ww(i+1,k) * column(QQ, k);
			//if(abs(weight) > eps)
		
		norm = inner_product(q2.begin(), q2.end(), q2.begin(), 0.0);
		norm = sqrt(norm);
		q2 /= norm;

		// Prepair for the next Lanczos step
		q0 = q1;
		q1 = q2;
	}
	if(i < nc)
		nEig = i;

	tridiag.resize(nEig, nEig + 1);
	for(i = 0; i < nEig; ++i) {
		for(k = 0; k < nEig; ++k)
			tridiag(i,k) = 0.0;
		tridiag(i,i) = al[i];
		if(i > 0)
			tridiag(i,i-1) = bt[i];
		if(i < nc - 1)
			tridiag(i,i+1) = bt[i+1];
	}
	tridiag.resize(nEig, nEig);

	std::vector<double> eig(nEig);
	matrix<double> eigVecTri(nc, nEig);
	matrix<double> eigVecSorted (nc, nEig);
	
	svdcmp(tridiag, eig, eigVecTri);
 
	std::vector<unsigned int> indx;
	sortWithIndex(eig, indx);

	std::vector<double> eigenval_withDuplicates(nEig); 
	for (i=0; i<nEig; i++) {
		eigenval_withDuplicates[i] = eig[indx[i]];
	}
	
	// Exclude duplicates
	double dist(1000.0);
	bool hasDuplicates(false);

	acc = 0.0;
	for(i=0; i<nEig; ++i) {
		if(i > 0) {
			if(eigenval_withDuplicates[i] < eps)
				dist = abs(eigenval_withDuplicates[i-1] - eigenval_withDuplicates[i]);
			else
				dist = abs(eigenval_withDuplicates[i-1]/ eigenval_withDuplicates[i] - 1.0);
		}
		if(dist > eps) {
			eigenval[icount] = eigenval_withDuplicates[i];
			acc += eigenval[icount];
			for(j=0; j<nc; j++)  
				eigVecSorted(j, icount) = eigVecTri(j,indx[i]);
			++icount;
		}
		else
			hasDuplicates = true;
	}

	eigenvec = prod(QQ, eigVecSorted);

	int nPC_out(min(nPC, icount));
	eigenvec_out.resize(eigenvec.size1(), nPC_out);
	eigenval_out.resize(nPC_out);
	for(i = 0; i < nPC_out; ++i) {
		column(eigenvec_out, i) = column(eigenvec, i);
		eigenval_out[i] = eigenval[i];
	}
	/*
	// Orthonormality check
	for(i = 0; i < nEig; ++i)
		for(j = 0; j < nEig; ++j) {
			vv = inner_product(column(eigVecSorted, j).begin(), column(eigVecSorted, j).end(), column(eigVecSorted, i).begin(), 0.0);
			if((i != j) && (abs(vv) > 10 *eps))
				eps = eps;
			if((i == j) && (abs(vv - 1) > 10 *eps))
				eps = eps;
		}
	// Orthonormality check
	for(i = 0; i < nEig; ++i)
		for(j = 0; j < nEig; ++j) {
			vv = inner_product(column(QQ, j).begin(), column(QQ, j).end(), column(QQ, i).begin(), 0.0);
			if((i != j) && (abs(vv) > 10 *eps))
				eps = eps;
			if((i == j) && (abs(vv - 1) > 10 *eps))
				eps = eps;
		}
	// Orthonormality check
	for(i = 0; i < nEig; ++i)
		for(j = 0; j < nEig; ++j) {
			vv = inner_product(column(eigenvec, j).begin(), column(eigenvec, j).end(), column(eigenvec, i).begin(), 0.0);
			if((i != j) && (abs(vv) > 10 *eps))
				eps = eps;
			if((i == j) && (abs(vv - 1) > 10 *eps))
				eps = eps;
		}
	*/
}
void svdcmp(matrix<double>& a, std::vector<double>& w, matrix<double>& v)
{
	int flag, i, its, j, jj, k, l, nm;
	int n(a.size2()), m(a.size1());
	double c,f,h,s,x,y,z;
	double anorm(0.0),g(0.0),scale(0.0);
	double volatile temp;   // See Saul Teukolsky's post on a convergence problem of this routine.
	std::vector<double> rv1(n+1);
	if (m < n)
		return;
		//nrerror("SVDCMP: You must augment A with extra zero rows");
	
	for (i=0; i<n; i++) {
		l =  i+ 1;
		rv1[i]=scale * g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a(k,i));
			if (scale) {
				for (k=i;k<m;k++) {
					a(k,i) /= scale;
					s += a(k,i)*a(k,i);
				}
				f=a(i, i);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a(i,i)=f-g;
				for (j=l; j<n; j++) {
					for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
					f=s/h;
					for (k=i;k<m;k++) a(k,j) += f*a(k,i);
				}
				for (k=i;k<m;k++) a(k,i) *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i < m && i != n-1) {
			for (k=l;k<n;k++) scale += fabs(a(i,k));
			if (scale) {
				for (k=l;k<n;k++) {
					a(i,k) /= scale;
					s += a(i,k)*a(i,k);
				}
				f=a(i,l);
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a(i,l)=f-g;
				for (k=l;k<n;k++) rv1[k]=a(i,k)/h;
				for (j=l;j<m;j++) {
					for (s=0.0,k=l;k<n;k++) s += a(j,k)*a(i,k);
					for (k=l;k<n;k++) a(j,k) += s*rv1[k];
				}
				for (k=l;k<n;k++) a(i,k) *= scale;
			}
		}
		anorm=max(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++)
					v(j,i)=(a(i,j)/a(i,l))/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a(i,k)*v(k,j);
					for (k=l;k<n;k++) v(k,j) += s*v(k,i);
				}
			}
			for (j=l;j<n;j++) v(i,j)=v(j,i)=0.0;
		}
		v(i,i)=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=min(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		if (i < n-1)
			for (j=l;j<n;j++) a(i,j)=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a(k,i)*a(k,j);
				f=(s/a(i,i))*g;
				for (k=i;k<m;k++) a(k,j) += f*a(k,i);
			}
			for (j=i;j<m;j++) a(j,i) *= g;
		} else {
			for (j=i;j<m;j++) a(j,i)=0.0;
		}
		++a(i,i);
	}
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=300;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				temp = fabs(rv1[l]) + anorm;  // See Saul Teukolsky's post on a convergence problem of this routine.
				if (temp == anorm) {
					flag=0;
					break;
				}
				temp = fabs(w[nm]) + anorm;  // See Saul Teukolsky's post on a convergence problem of this routine.
				if (temp == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					temp = fabs(f) + anorm;  // See Saul Teukolsky's post on a convergence problem of this routine.
					if ( temp == anorm ) break;
					g=w[i];
					h=PYTHAG(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s=(-f*h);
					for (j=0;j<m;j++) {
						y=a(j,nm);
						z=a(j,i);
						a(j,nm)=y*c+z*s;
						a(j,i)=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v(j,k)=(-v(j,k));
				}
				break;
			}
			//printf("Number of iterations in SVDCMP: %d\n",its);
			if (its == 300) 
				return;
				//nrerror("No convergence in 300 SVDCMP iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=PYTHAG(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=PYTHAG(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y=y*c;
				for (jj=0;jj<n;jj++) {
					x=v(jj,j);
					z=v(jj,i);
					v(jj,j)=x*c+z*s;
					v(jj,i)=z*c-x*s;
				}
				z=PYTHAG(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=(c*g)+(s*y);
				x=(c*y)-(s*g);
				for (jj=0;jj<m;jj++) {
					y=a(jj,j);
					z=a(jj,i);
					a(jj,j)=y*c+z*s;
					a(jj,i)=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	return;
}

void sortWithIndex(std::vector<double>& vec, std::vector<unsigned int>& index)
{
	unsigned int i, n(vec.size());
	std::vector< std::pair<int, double> > pairsVec;
	std::pair<int, double> indexValPair;

	for (i = 0; i < n; i++) {
		indexValPair.first = i;
		indexValPair.second = vec[i];
		pairsVec.push_back(indexValPair);
	}

	sort(pairsVec.begin(), pairsVec.end(), pairsValComp);
	index.resize(n);
	for (i = 0; i < n; i++) {
		index[i] = pairsVec[i].first;
	}
}
bool pairsValComp(const std::pair<int, double>& p1, const std::pair<int, double> p2)
{
	return p1.second > p2.second;
}

double orthog_check(matrix<double>& mat) {
// Check orthogonality of matrix columns. Returns average deviation from zero for off-diagonal elements
	double vv(0.0), mx(0.0);
	unsigned int i, j, nr(mat.size1()), nc(mat.size2());

	for(i = 0; i < nc; ++i) {
		matrix_column<matrix<double> > rw1 (mat, i);
		for(j = i; j < nc; ++j) {
			matrix_column<matrix<double> > rw2 (mat, j);
			vv = inner_product(rw1.begin(), rw1.end(), rw2.begin(), 0.0);
			if(i != j)
				mx += abs(vv);
		}
	}
	mx /= (nc * (nc + 1) / 2.0);
	return mx;
}
void gramm_shmidt_orth(matrix<double>& mat) { 
// Orthogonalize columns of mat
	double ss;
	unsigned int j, k, nc(mat.size2()), nr(mat.size1());
	std::vector<double> wts(nc - 1);

	for(j = 0; j < nc; ++j) {

		for(k = 0; k < j; ++k) 
			wts[k] = inner_product(column(mat, j).begin(), column(mat, j).end(), column(mat, k).begin(), 0.0) / nr;

		for(k = 0; k < j; ++k)	
			column(mat, j) -= wts[k] * column(mat, k);
		
		ss = inner_product(column(mat, j).begin(), column(mat, j).end(), column(mat, j).begin(), 0.0) / nr;
		ss = 1.0/sqrt(ss);
		
		column(mat, j) *= ss;
	}
}
namespace ublas = boost::numeric::ublas; 
bool InvertMatrix (const matrix<double>& input, matrix<double>& inverse) { 
    using namespace boost::numeric::ublas; 
    typedef permutation_matrix<std::size_t> pmatrix; 
    // create a working copy of the input 
    matrix<double> A(input); 
    // create a permutation matrix for the LU-factorization 
    pmatrix pm(A.size1()); 
    // perform LU-factorization 
    int res = lu_factorize(A,pm); 
    if( res != 0 )
        return false; 
    // create identity matrix of "inverse" 
    inverse.assign(ublas::identity_matrix<double>(A.size1())); 
    // backsubstitute to get the inverse 
    lu_substitute(A, pm, inverse); 
    return true; 
}
void calculate_weights(double decay, std::vector<date>& dts, std::vector<double>& wts, unsigned int sim_horizon) {
	// Input dates are ordered from past to present
	if(dts.size() < 2) {
		wts.resize(1);
		wts[0] = 1.0;
		return;
	}

	unsigned int n(dts.size() - sim_horizon);
	int i;
	double vv(0);
	wts.clear();
	wts.resize(n);

	if(decay < 1.0 && n > 1) {
		wts[n - 1] = 1.0;
		vv = 1.0;
		for(i = n-2; i >= 0; --i) {
			date_duration dt(dts[i + 1] - dts[i]);
			wts[i] = wts[i + 1] * pow(decay, dt.days());
			vv += wts[i];
		}
		vv = 1.0 / vv;

		for(i = 0; i < n; ++i) 
			wts[i] = sqrt(wts[i] * vv);

		vv = 0;
		for(i = 0; i < n; ++i)
			vv += wts[i] * wts[i];
		if(vv > 1)
			int tt = 0;
	}
	else {
		vv = sqrt(1.0 / n);
		for(i = 0; i < n; ++i)
			wts[i] = vv;
	}
}
double fit_normal_cdf_std(double x, double y) {
	// Finds standard deviation by fitting x and y to a function (normal cdf) by bisection
	// Should be generalized to a root search for any function Tanya

	double a(1.0), a1(0.001), a2(10), eps(1e-11), z1(1.0), z2(-1.0), z;
	unsigned int iter(0), maxIter(100);
	//normal norm_dist;
	//z1 =cdf(norm_dist, x / a1) - y;
	//z2 = cdf(norm_dist, x / a2) - y;
	z1 =cdf(x / a1) - y;
    z2 = cdf(x / a2) - y;

	while(iter < maxIter && abs(a1 - a2) > eps && z1 * z2 < 0) {
		
		a = 0.5 * (a1 + a2);
		//z = cdf(norm_dist, x / a) - y;
		z = cdf(x / a) - y;
		if(z * z1 < 0) {
			z2 = z;
			a2 = a;
		}
		else {
			z1 = z;
			a1 = a;
		}
		++iter;
	}

	return a;
}
int genIndepNorm(int m, std::vector<double>& vec) {
	
	boost::mt19937                     gener(1);
    boost::normal_distribution<double> normal_dist;   // Normal Distribution
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > rng(gener, normal_dist);
	rng.engine().seed(1);
	//rng.distribution().reset(); // to change the seed

	vec.resize(m);
    for (int i = 0; i < m; i++) 		
		vec[i] = rng();
	
	return 0;
}
double stdev(std::vector<double>& vec, std::vector<double>& wts) {
	double r(0.0), std(0.0), ave(0);
	unsigned int n(vec.size()), n1(n-1), i;
	if(n > 1) {
		for(i = 0; i < n; ++i)
			ave += vec[i];
		ave /= n;

		for(i = 0; i < n; ++i) 	{
			r = (vec[i] - ave) * wts[i];
			std += r * r;
		}
		std = sqrt(std);
	}
	return std;
}

bool regress(matrix<double>& X, matrix<double>& Y, matrix<double>& coef, std::vector<double>& vols_resid, std::vector<double>& R2) {
	// Each column of Y is regressed against X
	unsigned int i, j, n(X.size1()), m(X.size2()), k(Y.size2());
	std::vector<double> tt(m);

	if(X.size1() != Y.size1() || n == 0 || m == 0 || k == 0)
		return false;

	boost::numeric::ublas::vector<double> v(n), v1(n);
	matrix<double> aa(m, m), bb(m,k), invm(m,m), Y1(n, k), resid(n, k);
	double std(0);

	vols_resid.clear();
	R2.clear();
	coef.clear();
	vols_resid.resize(k);
	R2.resize(k);
	coef.resize(m, k);

	aa = prod(trans(X), X);
	for(i = 0; i < aa.size1(); ++i) {
		for(j = 0; j < aa.size2(); ++j)
			tt[j] = aa(i, j);
	}
	if(InvertMatrix(aa, invm)) {
		bb = prod(trans(X), Y);
		coef = prod(invm, bb);
		Y1 = prod(X, coef);
		resid = Y - Y1;
		for(i = 0; i < k; ++i) {
			vols_resid[i] = sqrt(inner_product(column(resid, i).begin(), column(resid, i).end(), column(resid, i).begin(), 0.0));
			std = sqrt(inner_product(column(Y, i).begin(), column(Y, i).end(), column(Y, i).begin(), 0.0));
			R2[i] = 1.0 - pow(vols_resid[i] / std, 2);
		}
	}
	return true;
}

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif 

double cdf(double x)
{
  double L, K, w ;
  
  double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
  double const a4 = -1.821255978, a5 = 1.330274429;

  L = fabs(x);
  K = 1.0 / (1.0 + 0.2316419 * L);
  w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

  if (x < 0 ){
    w= 1.0 - w;
  }
  return w;
}

void transform_2_empirical(std::vector<double>& empir_vec, std::vector<double>& norm_vec, std::vector<double>& res) {

	unsigned int i, j, nsim(norm_vec.size()), nrets(empir_vec.size());
	boost::numeric::ublas::vector<double> empir_prob(nrets), cdf_vec(nsim);
	double std(1.0), x1, x2, y1, y2, x, y, norm_mean, norm_std;
	
	normalize_vec(norm_vec, norm_mean, norm_std);
	sort(empir_vec.begin(), empir_vec.end());

	for(i = 0; i < nrets; ++i)
		empir_prob[i] = ( i + 1.0) / nrets;

	for(i = 0; i < nsim; ++i) {
		y = cdf(norm_vec[i]);
		j = 0;
		if(y <= empir_prob[j] || y >= empir_prob[nrets - 1]) // Simulated normal return is less than the smallest empirical return or greater than the largest empirical return - let it be
			res[i] = norm_vec[i];
		else {
			while(j < nrets && y > empir_prob[j] )
				j++;
			
			x1 = empir_vec[j - 1];
			y1 = empir_prob[j - 1];
			
			x2 = empir_vec[j];
			y2 = empir_prob[j];

			x = x1 +(x2 - x1) / (y2 - y1) * (y - y1);
			res[i] = x;
		}
	}
}
void normalize_vec(std::vector<double>& vec, double& mean, double& vol) {
	/*
	mean = std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
	for(int i = 0; i < vec.size(); ++i) 
		vec[i] -= mean;
	*/
	vol = sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0) / vec.size());
	for(int i = 0; i < vec.size(); ++i) 
		vec[i] /= vol;
}
bool interpolate_linear(std::vector<double>& x, std::vector<double>&y, std::vector<double>& arg, std::vector<double>& res) {
	std::map<double, double> mp;
	unsigned int i, n(x.size());

	if(x.size() != y.size())
		return false;

	for(i = 0; i < n; ++i) 
		mp.insert(std::make_pair(x[i], y[i]));

	interpolate_linear(mp, arg, res);
	return true;
}

void interpolate_linear(std::map<double,double>& mp, std::vector<double>& arg, std::vector<double>& res) {
	unsigned int i, m(arg.size());
	double val(0.0);
	std::map<double, double>::iterator it, jt;
	res.clear();
	res.resize(m, 0.0);

	if(mp.size() == 0)
		return;

	it = mp.begin();
	if(mp.size() == 1)
		for(i = 0; i < m; ++i) 
			res[i] = it->second;

	for(i = 0; i < m; ++i) {
		it = mp.upper_bound(arg[i]);
		if(it == mp.end())
			val = (--it)->second;
		else if(it == mp.begin())
			val = it->second;
		else {
			jt = it;
			--jt;
			val = jt->second + (it->second - jt->second) / (it->first - jt->first) * (arg[i] *  1.0 - jt->first);
		}
		res[i] = val;
	}
}
void interpolate_fwd_linear(std::map<double,double>& mp, std::vector<double>& arg, std::vector<double>& res) {
	unsigned int i, m(arg.size());
	double val(0.0), t1, t2, ff, y1, y2;
	std::map<double, double>::iterator it, jt;
	res.clear();
	res.resize(m, 0.0);

	if(mp.size() == 0)
		return;

	it = mp.begin();
	if(mp.size() == 1) {
		for(i = 0; i < m; ++i) 
			res[i] = it->second;
		return;
	}

	for(i = 0; i < m; ++i) {
		it = mp.upper_bound(arg[i]);
		if(it == mp.end())
			val = (--it)->second;
		else if(it == mp.begin())
			val = it->second;
		else {
			jt = it;
			--jt;
			t1 = jt->first;
			t2 = it->first;
			y1 = jt->second;
			y2 = it->second;
			ff = (y2 * t2 - y1 * t1) / (t2 - t1);
			val = ff + (y1 - y2) * t1 * t2 / (t2 - t1) / arg[i];
		}
		res[i] = val;
	}
}
void interpolate_quadratic(std::vector<double>& xi, std::vector<double>& yi, std::vector<double>& xo, std::vector<double>& yo) {
	unsigned int i, j, k, m, ni(xi.size()), no(xo.size()), N(4);
	double x0, y0, al, pp;
	std::vector<double> xx(N), yy(N);
	boost::numeric::ublas::vector<double> rh(N), res(N);
	matrix<double> aa(N, N), invm(N, N);

	yo.resize(no);

	if(yi.size() != ni)
		return;

	for(i = 0; i < no; i++) {
		x0 = xo[i];
		y0 = yo[i];
		j = 0;
		while(j < ni && x0 > xi[j])
			++j;
		
		if(j == 0 || j == ni)
			continue;

		if(j == 1) {
			al = -log(xi[1] / xi[0]) / log(yi[1]/yi[0]);
			y0 = yi[1] * pow(xi[1]/x0, 1.0 / al);
		}
		else if(j == ni - 1 ) {
			al = -log((1.0 - xi[ni - 2]) /(1.0 - xi[ni - 1])) / log(yi[ni - 2]/yi[ni - 1]) ;
			y0 = yi[ni - 2] * pow((1-xi[ni-2])/(1-x0), 1.0 / al);
		}
		else  {
			for(k = 0; k < N; ++k) {
				rh[k] = yi[j - 2 + k]; 
				xx[k] = xi[j - 2 + k];
				aa(k, 0) = 1.0;
			}
				
			for(k = 1; k < N; ++k) {
				for(m = 0; m < N; ++m)
					aa(m, k) = aa(m, k - 1) * xx[m];
			}
			if(InvertMatrix(aa, invm)) {
				res = prod(invm, rh);
			}
			y0 = res[0];
			pp = 1.0;
			for(k = 1; k < N; ++k) {
				pp *= x0;
				y0 += pp * res[k];
			}
		}
		yo[i] = y0;
	}
}
/*
void transform_2_empirical(matrix_column<matrix<double> > empir_vec, matrix_column<matrix<double> > norm_vec, matrix_column<matrix<double> >& res) {

	unsigned int i, j, nsim(norm_vec.size()), nrets(empir_vec.size());
	boost::numeric::ublas::vector<double> empir_prob(nrets), cdf_vec(nsim);
	double std(1.0), x1, x2, y1, y2, x, y, mean;

	mean = std::accumulate(norm_vec.begin(), norm_vec.end(), 0.0) / norm_vec.size();
	std = sqrt(inner_product(norm_vec.begin(), norm_vec.end(), norm_vec.begin(), 0.0) / nsim - mean * mean);
	norm_vec /= std;

	sort(empir_vec.begin(), empir_vec.end());
	for(i = 0; i < nrets; ++i)
		empir_prob[i] = ( i + 1.0) / nrets;

	for(i = 0; i < nsim; ++i) {
		y = cdf(norm_vec[i]);
		j = 0;
		if(y <= empir_prob[j] || y >= empir_prob[nrets - 1]) // Simulated normal return is less than the smallest empirical return or greater than the largest empirical return - let it be
			res[i] = norm_vec[i];
		else {
			while(j < nrets && y > empir_prob[j] )
				j++;
			
			x1 = empir_vec[j - 1];
			y1 = empir_prob[j - 1];
			
			x2 = empir_vec[j];
			y2 = empir_prob[j];
			x = x1 +(x2 - x1) / (y2 - y1) * (y - y1);
			res[i] = x;
		}
	}
}
*/
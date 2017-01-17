#ifndef CCARSHIFT_H
#define CCARSHIFT_H

class CCAR_shift {

public:
	CCAR_shift()  : _isLevel(true), _isMult(true), _s_name("") {
		_terms.resize(1);
		_terms[0] = 0.0;
		_shifts.resize(1,1);
		_shifts(0, 0) = 1.0;
		};
	CCAR_shift(const CCAR_shift& vv) : _isLevel(vv._isLevel), _isMult(vv._isMult), _s_name(vv._s_name), _terms(vv._terms), _shifts(vv._shifts) {};
	~CCAR_shift() {};
	void setName(std::string nm) {_s_name = nm;};
	void setIsMult(bool vv) { _isMult = vv;};
	void setIsLevel(bool vv) { _isLevel = vv;};
	void setTerms(const std::vector<double>& vec) {_terms = vec;};
	void setShifts(matrix<double>& sh) {_shifts = sh;};
	std::string getName() const {return _s_name;};
	bool					getIsLevel() {return _isLevel;};
	bool					getIsMult() {return _isMult;};
	const std::vector<double>&	getTerms() const {return _terms;};
	const matrix<double>&			getShifts() const {return _shifts;};

private:
	bool				_isLevel; // to shift lecels or vols
	bool				_isMult; // it multiplicative or additive
	std::string			_s_name;
	std::vector<double> _terms;
	matrix<double>		_shifts;
};

#endif
#ifndef SHIFTSCENINFO_H
#define SHIFTSCENINFO_H

class ShiftScenInfo {
public:
	ShiftScenInfo(std::string nm1, std::string nm2, unsigned int n1, unsigned int n2, ShiftType sh, double v);
	ShiftScenInfo() : _s_name1(""), _s_name2(""), _tenor1(0), _tenor2(0), _shiftType(Base), _shift(0.0) {};
	~ShiftScenInfo() {};
	std::string		getName1() {return _s_name1;};
	std::string		getName2() {return _s_name2;};
	unsigned int	getTenor1() {return _tenor1;};
	unsigned int	getTenor2() {return _tenor2;};
	ShiftType		getShiftType() {return _shiftType;};
	double			getShift() {return _shift;};

	void setName1(std::string nm) {_s_name1 = nm;};
	void setName2(std::string nm) {_s_name2 = nm;};
	void setTenor1(unsigned int n) {_tenor1 = n;};
	void setTenor2(unsigned int n) {_tenor2 = n;};
	void setShiftType(ShiftType type) {_shiftType = type;};
	void setShift(double sh) {_shift = sh;};

private:
	std::string _s_name1;
	std::string	_s_name2;
	unsigned int _tenor1;
	unsigned int _tenor2;
	ShiftType	 _shiftType;
	double		 _shift;
};
#endif SHIFTSCENINFO_H
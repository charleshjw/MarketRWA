#ifndef COREFACTOR_H
#define COREFACTOR_H

class CoreFactor : public Factor {

public:
	CoreFactor(const Factor& ff);
	~CoreFactor() {};

private:
	bool _isMultiplicative;
	std::vector<double> tenors;
	matrix<double> shifts;
};

#endif COREFACTOR_H
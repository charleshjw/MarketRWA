#include "stdafx.h"
//#include <fstream>

using namespace std;

ErrorLog::ErrorLog(string file_name) :_file_name(file_name) {
	// Filename should be err_dir + "err_"
	if(file_name !="")
		_out_file.open(file_name.c_str());
}
ErrorLog::~ErrorLog() {
	if(!_out_file.fail())
		_out_file.close();
}
void ErrorLog:: print() {
	if(!_out_file.fail()) {
		unsigned int i(0);
		for(i = 0; i < _errors.size(); ++i) {
			_out_file << _errors[i] << endl;
		}
	}
}
void ErrorLog::add_line(string line) {
	if(line != "")
		_errors.push_back(line);
}
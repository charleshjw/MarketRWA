#ifndef ERRORLOG_H
#define ERRORLOG_H

class ErrorLog {
public:
	ErrorLog(std::string file_name);
	~ErrorLog();
	void add_line(std::string line);
	void print();
private:
	std::vector<std::string> _errors;
	std::string _file_name;
	std::ofstream _out_file;
};
#endif ERRORLOG
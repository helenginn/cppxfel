#include <string>
#include <stdio.h>
#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <vector>

std::string i_to_str(int val);
std::string f_to_str(double val);

bool replace(std::string& str, const std::string& from, const std::string& to);

std::string getBaseFilename(std::string filename);
std::string getFilename(std::string filename);
std::string getPath(std::string filename);

std::vector<std::string> glob(std::string globString);

unsigned long factorial(unsigned long n);
unsigned int choose(unsigned long n, unsigned long choose);
double proportion(int n);

// trim from start
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}
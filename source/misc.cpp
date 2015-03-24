#include <string>
#include <sstream>
#include "misc.h"
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>

string f_to_str(double val)
{
	ostringstream ss;
	ss << val;
	string temp = ss.str();

	return temp;
}

string i_to_str(int val)
{
	ostringstream ss;
	ss << val;
	string temp = ss.str();
	
	return temp;
}

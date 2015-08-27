#include <string>
#include <sstream>
#include "misc.h"
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>

std::string f_to_str(double val)
{
	std::ostringstream ss;
	ss << val;
	std::string temp = ss.str();

	return temp;
}

std::string i_to_str(int val)
{
	std::ostringstream ss;
	ss << val;
	std::string temp = ss.str();
	
	return temp;
}

bool replace(std::string& str, const std::string& from, const std::string& to)
{
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::string getBaseFilename(std::string filename)
{
    std::string fName = getFilename(filename);
    size_t pos = fName.rfind(".");
    if(pos == std::string::npos)  //No extension.
        return fName;
    
    if(pos == 0)    //. is at the front. Not an extension.
        return fName;
    
    return fName.substr(0, pos);
}

std::string getFilename(std::string filename)
{
    size_t pos = filename.rfind("/");
    if(pos == std::string::npos)  //No path.
        return filename;
    
    return filename.substr(pos + 1, filename.length());
}

std::string getPath(std::string filename)
{
    size_t pos = filename.rfind("/");
    if(pos == std::string::npos)  //No path.
        return filename;
    
    return filename.substr(0, pos + 1);
}
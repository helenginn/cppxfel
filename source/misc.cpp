#include <string>
#include <sstream>
#include "misc.h"
#include <iostream>
#include <sys/wait.h>
#include <unistd.h>
#include <math.h>
#include <glob.h>
#include <vector>


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

unsigned long factorial(unsigned long n)
{
    unsigned long ret = 1;
    for(unsigned int i = 1; i <= n; ++i)
        ret *= i;
    return ret;
}

unsigned int choose(unsigned long n, unsigned long choose)
{
    if (n > 18) n = 18;
    
    unsigned long nFac = factorial(n);
    unsigned long cFac = factorial(choose);
    unsigned long ncFac = factorial(n - choose);
    
    unsigned int value = (unsigned int)(nFac / (cFac * ncFac));
    
    return value;
}

double proportion(int n)
{
    double prop = n / pow(2, n);
    
    return prop;
}


std::vector<std::string> glob(std::string globString)
{
    std::vector<std::string> globs;
    
    glob_t globbuf;
    int err = glob(globString.c_str(), 0, NULL, &globbuf);
    if(err == 0)
    {
        for (size_t i = 0; i < globbuf.gl_pathc; i++)
        {
            globs.push_back(std::string(globbuf.gl_pathv[i]));
        }
        
        globfree(&globbuf);
    }
    else
    {
        switch (err)
        {
            case GLOB_NOMATCH:
                std::cout << "GlOB_NOMATCH" << std::endl;
                break;
            case GLOB_ABORTED:
                std::cout << "GLOB_ABORTED" << std::endl;
                break;
            case GLOB_NOSPACE:
                std::cout << "GLOB_NOSPACE" << std::endl;
                break;

            default:
                break;
        }
        
    }
    
    return globs;
}

void trim(std::string& str)
{
    std::string::size_type pos = str.find_last_not_of(' ');
    if(pos != std::string::npos) {
        str.erase(pos + 1);
        pos = str.find_first_not_of(' ');
        if(pos != std::string::npos) str.erase(0, pos);
            }
    else str.erase(str.begin(), str.end());
        }
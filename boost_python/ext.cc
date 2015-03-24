/*
 * flex_ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <iostream>
#include <string>
#include "../source/PythonExt.cpp"

namespace cppxfel { namespace boost_python {

	using namespace boost::python;

	void printString(std::string stringToPrint)
	{
		std::cout << stringToPrint << std::endl;
	}

	void cppxfelScript(std::string scriptFile)
	{
		runScriptFromPython(scriptFile);
	}

	BOOST_PYTHON_MODULE(cppxfel_ext)
	{
 		def ("printString", &printString);
 		def ("cppxfel", &cppxfelScript);
 	}
}

}
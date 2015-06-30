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
#include "../source/MtzRefiner.h"
#include "../source/parameters.h"
#include <iostream>
#include <string>
#include "../source/PythonExt.h"
#include "../source/InputFileParser.h"
#include <boost/container/vector.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "../source/MtzManager.h"

namespace cppxfel { namespace boost_python {

	using namespace boost::python;

	template<class T>
	boost::python::list vector_to_py_list(const vector<T>& v)
	{
		boost::python::object get_iter = boost::python::iterator<vector<T> >();
		boost::python::object iter = get_iter(v);
		boost::python::list l(iter);
		return l;
	}
	
	boost::python::list getMtzs(InputFileParser parser)
	{
		vector<MtzPtr> mtzs = parser.getRefiner()->getMtzManagers();
		boost::python::list list = vector_to_py_list(mtzs);
		return list;
	}

	void cppxfelScript(std::string scriptFile)
	{
		runScriptFromPython(scriptFile);
	}

	BOOST_PYTHON_MODULE(cppxfel_ext)
	{
 		def ("run", &cppxfelScript);
 		
 		class_<MtzManager, MtzPtr, boost::noncopyable>("MtzManager", no_init)
 			.def("gridSearch", &MtzManager::gridSearch)
 			.def("refCorrelation", &MtzManager::getRefCorrelation)
 			.def("applyUnrefinedPartiality", &MtzManager::applyUnrefinedPartiality)
 		;
 		
 		class_<std::vector<MtzPtr> >("MtzArray")
			.def(vector_indexing_suite<vector<MtzPtr> >()); 		
 		
		class_<InputFileParser>("Parser", init<std::string>())
			.def("parse", &InputFileParser::parse)
			.def("refine", &InputFileParser::refine)
			.def("mtzs", &getMtzs)
		;
 	}
}

}
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
#include <dials/model/data/shoebox.h>
#include <dials/algorithms/centroid/centroid.h>

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

        template<class T>
        std::vector<T> py_list_to_vector(boost::python::list theList)
        {
                std::vector<T> vec;

                for (int i = 0; i < len(theList); ++i)
    {
        vec.push_back(boost::python::extract<T>(theList[i]));
    }

                return vec;
        }

        void cppxfelScript(std::string scriptFile)
        {
                runScriptFromPython(scriptFile);
        }

        BOOST_PYTHON_MODULE(cppxfel_ext)
        {
                def ("run", &cppxfelScript);
                def ("runCommandLineArgs", &runCommandLine);
        }
}

}

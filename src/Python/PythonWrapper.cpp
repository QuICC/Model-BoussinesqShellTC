/** 
 * @file PythonWrapper.cpp
 * @brief Source of Python interpreter wrapper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Python/PythonWrapper.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

   PyObject* PythonWrapper::mpModule = NULL;

   PyObject* PythonWrapper::mpFunc = NULL;

   void PythonWrapper::init()
   {
      // Initialize the python interpreter
      Py_Initialize();

      // Setup the search path
      PyObject* sysPath = PySys_GetObject((char*)"path");
      PyList_Append(sysPath, PyUnicode_FromString("."));
   }

   void PythonWrapper::import(const std::string& module)
   {
      // Get string object for module name
      PyObject* pName;
      pName = PyUnicode_FromString(module.c_str());

      // Import module
      PythonWrapper::mpModule = PyImport_Import(pName);

      // Release pName
      Py_DECREF(pName);

      if(PythonWrapper::mpModule == NULL)
      {
         PyErr_Print();
         throw Exception("Python module import error!");
      }
   }

   void PythonWrapper::setFunction(const std::string& func)
   {
      // Get Python function object
      PythonWrapper::mpFunc = PyObject_GetAttrString(PythonWrapper::mpModule, func.c_str());

      // Check for successfully loading function
      if(! PythonWrapper::mpFunc || PyCallable_Check(PythonWrapper::mpFunc))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw Exception("Python function loading error!");
      }
   }

   void PythonWrapper::finalize()
   {
      // Clean up
      if(PythonWrapper::mpFunc != NULL)
      {
         Py_DECREF(PythonWrapper::mpFunc);
      }
      if(PythonWrapper::mpModule != NULL)
      {
         Py_DECREF(PythonWrapper::mpModule);
      }

      // Finalize the python interpreter
      Py_Finalize();
   }

   PythonWrapper::PythonWrapper()
   {
   }

   PythonWrapper::~PythonWrapper()
   {
   }

}

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
#include "Python/PythonConfig.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/HumanToId.hpp"

namespace GeoMHDiSCC {

   PyObject* PythonWrapper::mpModule = NULL;

   PyObject* PythonWrapper::mpFunc = NULL;

   void PythonWrapper::init()
   {
      // Initialize the python interpreter
      Py_Initialize();

      // Setup the search path
      PyObject* sysPath = PySys_GetObject((char*)"path");
      PyList_Append(sysPath, PyUnicode_FromString(GEOMHDISCC_PYTHON_DIR));
   }

   void PythonWrapper::import(const std::string& module)
   {
      // Cleanup
      if(PythonWrapper::mpModule != NULL)
      {
         Py_CLEAR(PythonWrapper::mpModule);
      }

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
      // Cleanup
      if(PythonWrapper::mpFunc != NULL)
      {
         Py_CLEAR(PythonWrapper::mpFunc);
      }

      // Get Python function object
      PythonWrapper::mpFunc = PyObject_GetAttrString(PythonWrapper::mpModule, func.c_str());

      // Check for successfully loading function
      if(! (PythonWrapper::mpFunc && PyCallable_Check(PythonWrapper::mpFunc)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw Exception("Python function loading error!");
      }
   }

   PyObject* PythonWrapper::callFunction()
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(PythonWrapper::mpFunc, NULL);

      return pValue;
   }

   PyObject* PythonWrapper::callFunction(PyObject* pArgs)
   {
      PyObject *pValue;

      pValue = PyObject_CallObject(PythonWrapper::mpFunc, pArgs);

      return pValue;
   }

   PyObject* PythonWrapper::makeTuple(const ArrayI& val)
   {
      PyObject *pTuple, *pValue;

      pTuple = PyTuple_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = PyLong_FromLong(val(i));
         PyTuple_SetItem(pTuple, i, pValue); 
      }

      return pTuple;
   }

   PyObject* PythonWrapper::makeTuple(const std::vector<MHDFloat>& val)
   {
      PyObject *pTuple, *pValue;

      pTuple = PyTuple_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = PyFloat_FromDouble(val.at(i));
         PyTuple_SetItem(pTuple, i, pValue); 
      }

      return pTuple;
   }

   PyObject* PythonWrapper::makeDict(const std::vector<std::string>& key, const std::vector<MHDFloat>& val)
   {
      PyObject *pDict, *pKey, *pValue;

      pDict = PyDict_New();
      for(unsigned int i = 0; i < key.size(); i++)
      {
         pKey = PyUnicode_FromString(key.at(i).c_str());
         pValue = PyFloat_FromDouble(val.at(i));
         PyDict_SetItem(pDict, pKey, pValue);
      }

      return pDict;
   }

   PyObject* PythonWrapper::makeDict(const std::map<std::string, int>& map)
   {
      PyObject *pDict, *pKey, *pValue;

      pDict = PyDict_New();
      std::map<std::string,int>::const_iterator mapIt;
      for(mapIt = map.begin(); mapIt != map.end(); ++mapIt)
      {
         pKey = PyUnicode_FromString(mapIt->first.c_str());
         pValue = PyLong_FromLong(mapIt->second);
         PyDict_SetItem(pDict, pKey, pValue);
      }

      return pDict;
   }

   void PythonWrapper::getList(std::vector<bool> &rList, PyObject *pList)
   {
      PyObject *pValue;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         bool isPeriodic = PyObject_IsTrue(pValue);

         rList.push_back(isPeriodic);
      }
   }

   void PythonWrapper::getList(std::vector<NonDimensional::Id> &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyUnicode_AsASCIIString(pValue);
         NonDimensional::Id nd = IoTools::HumanToId::toNd(std::string(PyBytes_AsString(pTmp)));

         rList.push_back(nd);
      }
   }

   void PythonWrapper::getList(std::vector<PhysicalNames::Id> &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyUnicode_AsASCIIString(pValue);
         PhysicalNames::Id phys = IoTools::HumanToId::toPhys(std::string(PyBytes_AsString(pTmp)));

         rList.push_back(phys);
      }
   }

   void PythonWrapper::getList(std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> > &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp, *pTmp2;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyTuple_GetItem(pValue,0);
         pTmp2 = PyUnicode_AsASCIIString(pTmp);
         PhysicalNames::Id phys = IoTools::HumanToId::toPhys(std::string(PyBytes_AsString(pTmp2)));

         pTmp = PyTuple_GetItem(pValue,1);
         pTmp2 = PyUnicode_AsASCIIString(pTmp);
         FieldComponents::Spectral::Id comp = IoTools::HumanToId::toComp(std::string(PyBytes_AsString(pTmp2)));

         rList.push_back(std::make_pair(phys, comp));
      }
   }

   void PythonWrapper::fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat) 
   {
      PyObject *pArgs, *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));

      // Convert Python matrix into triplets
      pArgs = PyTuple_New(1);
      PyTuple_SetItem(pArgs, 0, pPyMat);
      PythonWrapper::setFunction((char *)"triplets");
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      long int len = PyList_Size(pValue);
      std::vector<Triplet> triplets;
      triplets.reserve(len);
      long int row;
      long int col;
      MHDFloat val;
      for(int i = 0; i < len; i++)
      {
         pTmp = PyList_GetItem(pValue, i);
         row = PyLong_AsLong(PyTuple_GetItem(pTmp,0));
         col = PyLong_AsLong(PyTuple_GetItem(pTmp,1));
         val = PyFloat_AsDouble(PyTuple_GetItem(pTmp,2));
         triplets.push_back(Triplet(row,col,val));
      }
      Py_DECREF(pValue);

      // Build matrix
      rMatrix.resize(rows,cols);
      rMatrix.setFromTriplets(triplets.begin(), triplets.end());
   }

   void PythonWrapper::fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat) 
   {
      PyObject *pArgs, *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));

      // Convert Python matrix into triplets
      pArgs = PyTuple_New(1);
      PyTuple_SetItem(pArgs, 0, pPyMat);
      PythonWrapper::setFunction((char *)"triplets");
      pValue = PythonWrapper::callFunction(pArgs);
      Py_DECREF(pArgs);

      long int len = PyList_Size(pValue);
      std::vector<Triplet> realTriplets;
      std::vector<Triplet> imagTriplets;
      realTriplets.reserve(len);
      imagTriplets.reserve(len);
      long int row;
      long int col;
      MHDFloat val;
      for(int i = 0; i < len; i++)
      {
         pTmp = PyList_GetItem(pValue, i);
         row = PyLong_AsLong(PyTuple_GetItem(pTmp,0));
         col = PyLong_AsLong(PyTuple_GetItem(pTmp,1));
         if(PyComplex_Check(PyTuple_GetItem(pTmp,2)))
         {
            val = PyComplex_RealAsDouble(PyTuple_GetItem(pTmp,2));
            realTriplets.push_back(Triplet(row,col,val));
            val = PyComplex_ImagAsDouble(PyTuple_GetItem(pTmp,2));
            imagTriplets.push_back(Triplet(row,col,val));
         } else
         {
            val = PyFloat_AsDouble(PyTuple_GetItem(pTmp,2));
            realTriplets.push_back(Triplet(row,col,val));
         }
      }
      Py_DECREF(pValue);

      // Build matrix
      rMatrix.real().resize(rows,cols);
      rMatrix.real().setFromTriplets(realTriplets.begin(), realTriplets.end());

      if(imagTriplets.size() > 0)
      {
         rMatrix.imag().resize(rows,cols);
         rMatrix.imag().setFromTriplets(imagTriplets.begin(), imagTriplets.end());
      } else
      {
         rMatrix.imag().resize(0,0);
      }
   }

   void PythonWrapper::cleanup()
   {
      // Clean up
      if(PythonWrapper::mpFunc != NULL)
      {
         Py_CLEAR(PythonWrapper::mpFunc);
      }
      if(PythonWrapper::mpModule != NULL)
      {
         Py_CLEAR(PythonWrapper::mpModule);
      }
   }

   void PythonWrapper::finalize()
   {
      // Clean up
      if(PythonWrapper::mpFunc != NULL)
      {
         Py_CLEAR(PythonWrapper::mpFunc);
      }
      if(PythonWrapper::mpModule != NULL)
      {
         Py_CLEAR(PythonWrapper::mpModule);
      }

      // Finalize
      Py_Finalize();
   }

   PythonWrapper::PythonWrapper()
   {
   }

   PythonWrapper::~PythonWrapper()
   {
   }

}

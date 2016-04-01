/** 
 * @file PythonTools.cpp
 * @brief Source of the tools to work with Python objects
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
#include "Python/PythonTools.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoTools/HumanToId.hpp"

namespace GeoMHDiSCC {

   PyObject* PythonTools::makeTuple(const ArrayI& val)
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

   PyObject* PythonTools::makeTuple(const std::vector<MHDFloat>& val)
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

   PyObject* PythonTools::makeList(const std::vector<int>& val)
   {
      PyObject *pList, *pValue;

      pList = PyList_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = PyLong_FromLong(val.at(i));
         PyList_SetItem(pList, i, pValue); 
      }

      return pList;
   }

   PyObject* PythonTools::makeList(const std::vector<std::vector<int> >& val)
   {
      PyObject *pList, *pValue;

      pList = PyList_New(val.size());
      for(unsigned int i = 0; i < val.size(); i++)
      {
         pValue = PythonTools::makeList(val.at(i));
         PyList_SetItem(pList, i, pValue); 
      }

      return pList;
   }

   PyObject* PythonTools::makeDict(const std::vector<std::string>& key, const std::vector<MHDFloat>& val)
   {
      PyObject *pDict, *pKey, *pValue;

      pDict = PyDict_New();
      for(unsigned int i = 0; i < key.size(); i++)
      {
         pKey = PyUnicode_FromString(key.at(i).c_str());
         pValue = PyFloat_FromDouble(val.at(i));
         PyDict_SetItem(pDict, pKey, pValue);
         Py_DECREF(pValue);
         Py_DECREF(pKey);
      }

      return pDict;
   }

   PyObject* PythonTools::makeDict(const std::map<std::string, int>& map)
   {
      PyObject *pDict, *pKey, *pValue;

      pDict = PyDict_New();
      std::map<std::string,int>::const_iterator mapIt;
      for(mapIt = map.begin(); mapIt != map.end(); ++mapIt)
      {
         pKey = PyUnicode_FromString(mapIt->first.c_str());
         pValue = PyLong_FromLong(mapIt->second);
         PyDict_SetItem(pDict, pKey, pValue);
         Py_DECREF(pValue);
         Py_DECREF(pKey);
      }

      return pDict;
   }

   void PythonTools::getList(std::vector<bool> &rList, PyObject *pList)
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

   void PythonTools::getList(std::vector<NonDimensional::Id> &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyUnicode_AsASCIIString(pValue);
         NonDimensional::Id nd = IoTools::HumanToId::toNd(std::string(PyBytes_AsString(pTmp)));
         Py_DECREF(pTmp);

         rList.push_back(nd);
      }
   }

   void PythonTools::getList(std::vector<PhysicalNames::Id> &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyUnicode_AsASCIIString(pValue);
         PhysicalNames::Id phys = IoTools::HumanToId::toPhys(std::string(PyBytes_AsString(pTmp)));
         Py_DECREF(pTmp);

         rList.push_back(phys);
      }
   }

   void PythonTools::getList(std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> > &rList, PyObject *pList)
   {
      PyObject *pValue, *pTmp, *pTmp2;

      long int len = PyList_Size(pList);
      for(int i = 0; i < len; ++i)
      {
         pValue = PyList_GetItem(pList, i);

         pTmp = PyTuple_GetItem(pValue,0);
         pTmp2 = PyUnicode_AsASCIIString(pTmp);
         PhysicalNames::Id phys = IoTools::HumanToId::toPhys(std::string(PyBytes_AsString(pTmp2)));
         Py_DECREF(pTmp2);

         pTmp = PyTuple_GetItem(pValue,1);
         pTmp2 = PyUnicode_AsASCIIString(pTmp);
         FieldComponents::Spectral::Id comp = IoTools::HumanToId::toComp(std::string(PyBytes_AsString(pTmp2)));
         Py_DECREF(pTmp2);

         rList.push_back(std::make_pair(phys, comp));
      }
   }

   PythonTools::PythonTools()
   {
   }

   PythonTools::~PythonTools()
   {
   }

}

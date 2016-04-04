/** 
 * @file PythonModelWrapper.cpp
 * @brief Source of the model Python interpreter wrapper
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
#include "Python/PythonModelWrapper.hpp"

// Project includes
//
#include "Python/PythonConfig.hpp"
#include "Python/PythonCoordinator.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/HumanToId.hpp"

namespace GeoMHDiSCC {

   PyObject* PythonModelWrapper::mpModule = NULL;

   PyObject* PythonModelWrapper::mpFunc = NULL;

   PyObject* PythonModelWrapper::mpModel = NULL;

   PyObject* PythonModelWrapper::mpMethod = NULL;

   void PythonModelWrapper::init()
   {
      PythonCoordinator::init();
   }

   void PythonModelWrapper::import(const std::string& module)
   {
      // Cleanup
      if(PythonModelWrapper::mpModule != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpModule);
      }

      // Get string object for module name
      PyObject* pName;
      pName = PyUnicode_FromString(("geomhdiscc.model."+module).c_str());

      // Import module
      PythonModelWrapper::mpModule = PyImport_Import(pName);

      // Release pName
      Py_DECREF(pName);

      if(PythonModelWrapper::mpModule == NULL)
      {
         PyErr_Print();
         throw Exception("Python module import error!");
      }
   }

   void PythonModelWrapper::createModel(const std::string& model)
   {
      // Cleanup
      if(PythonModelWrapper::mpModel != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpModel);
      }

      // Create model object

      PythonModelWrapper::setFunction(model);
      PythonModelWrapper::mpModel = PythonModelWrapper::callFunction();

      if(PythonModelWrapper::mpModel == NULL)
      {
         PyErr_Print();
         throw Exception("Python model creation error!");
      }

      #ifdef GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
      // Set the galerkin flag for the python dispatchers
      PyObject_SetAttrString(PythonModelWrapper::mpModel, (char *)"use_galerkin", Py_True);
      #endif //GEOMHDISCC_BOUNDARYMETHOD_GALERKIN
   }

   void PythonModelWrapper::setFunction(const std::string& func)
   {
      // Cleanup
      if(PythonModelWrapper::mpFunc != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpFunc);
      }

      // Get Python function object
      PythonModelWrapper::mpFunc = PyObject_GetAttrString(PythonModelWrapper::mpModule, func.c_str());

      // Check for successfully loading function
      if(! (PythonModelWrapper::mpFunc && PyCallable_Check(PythonModelWrapper::mpFunc)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw Exception("Python function loading error!");
      }
   }

   void PythonModelWrapper::setFunction(const std::string& func, const std::string& submod)
   {
      // Cleanup
      if(PythonModelWrapper::mpFunc != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpFunc);
      }

      // Get Python function object
      PyObject *pTmp = PyObject_GetAttrString(PythonModelWrapper::mpModule, submod.c_str());
      if(pTmp)
      {
         PythonModelWrapper::mpFunc = PyObject_GetAttrString(pTmp, func.c_str());
         Py_DECREF(pTmp);

         // Check for successfully loading function
         if(! (PythonModelWrapper::mpFunc && PyCallable_Check(PythonModelWrapper::mpFunc)))
         {
            if(PyErr_Occurred())
            {
               PyErr_Print();
            }
            throw Exception("Python function loading error!");
         }
      } else
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
      }
   }

   void PythonModelWrapper::setMethod(const std::string& method)
   {
      // Cleanup
      if(PythonModelWrapper::mpMethod != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpMethod);
      }

      // Get model method object
      PythonModelWrapper::mpMethod = PyObject_GetAttrString(PythonModelWrapper::mpModel, method.c_str());

      // Check for successfully loaded method
      if(! (PythonModelWrapper::mpMethod && PyCallable_Check(PythonModelWrapper::mpMethod)))
      {
         if(PyErr_Occurred())
         {
            PyErr_Print();
         }
         throw Exception("Python model method loading error!");
      }
   }

   PyObject* PythonModelWrapper::callFunction()
   {
      PyObject *pValue;

      pValue = PyEval_CallObject(PythonModelWrapper::mpFunc, NULL);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw Exception("Python function call error!");
      }

      return pValue;
   }

   PyObject* PythonModelWrapper::callFunction(PyObject* pArgs)
   {
      PyObject *pValue;

      pValue = PyEval_CallObject(PythonModelWrapper::mpFunc, pArgs);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw Exception("Python function with arguments call error!");
      }

      return pValue;
   }

   PyObject* PythonModelWrapper::callMethod()
   {
      PyObject *pValue;

      pValue = PyEval_CallObject(PythonModelWrapper::mpMethod, NULL);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw Exception("Python method call error!");
      }

      return pValue;
   }

   PyObject* PythonModelWrapper::callMethod(PyObject* pArgs)
   {
      PyObject *pValue;

      pValue = PyEval_CallObject(PythonModelWrapper::mpMethod, pArgs);

      if(PyErr_Occurred())
      {
         PyErr_Print();
         throw Exception("Python method with arguments call error!");
      }

      return pValue;
   }

   void PythonModelWrapper::fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat) 
   {
      PyObject *pArgs, *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Convert Python matrix into triplets
      pArgs = PyTuple_New(1);
      Py_INCREF(pPyMat);
      PyTuple_SetItem(pArgs, 0, pPyMat);
      PythonModelWrapper::setFunction("triplets", "utils");
      pValue = PythonModelWrapper::callFunction(pArgs);
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

   void PythonModelWrapper::fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat) 
   {
      PyObject *pArgs, *pValue, *pTmp;

      // Get matrix size
      pValue = PyObject_GetAttrString(pPyMat, (char *)"shape");
      long int rows = PyLong_AsLong(PyTuple_GetItem(pValue, 0));
      long int cols = PyLong_AsLong(PyTuple_GetItem(pValue, 1));
      Py_DECREF(pValue);

      // Convert Python matrix into triplets
      pArgs = PyTuple_New(1);
      Py_INCREF(pPyMat);
      PyTuple_SetItem(pArgs, 0, pPyMat);
      PythonModelWrapper::setFunction("triplets", "utils");
      pValue = PythonModelWrapper::callFunction(pArgs);
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
      rMatrix.imag().resize(rows,cols);
      rMatrix.real().setFromTriplets(realTriplets.begin(), realTriplets.end());

      if(imagTriplets.size() > 0)
      {
         rMatrix.imag().setFromTriplets(imagTriplets.begin(), imagTriplets.end());
      }
   }

   void PythonModelWrapper::cleanup()
   {
      // Clean up
      if(PythonModelWrapper::mpFunc != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpFunc);
      }
      if(PythonModelWrapper::mpMethod != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpMethod);
      }
   }

   void PythonModelWrapper::finalize()
   {
      // Clean up
      PythonModelWrapper::cleanup();

      // Clear model and module
      if(PythonModelWrapper::mpModel != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpModel);
      }
      if(PythonModelWrapper::mpModule != NULL)
      {
         Py_CLEAR(PythonModelWrapper::mpModule);
      }

      // Finalize
      PythonCoordinator::finalize();
   }

   PythonModelWrapper::PythonModelWrapper()
   {
   }

   PythonModelWrapper::~PythonModelWrapper()
   {
   }

}

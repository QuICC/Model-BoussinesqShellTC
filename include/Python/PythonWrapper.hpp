/**
 * @file PythonWrapper.hpp
 * @brief Small wrapper for python embedding
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PYTHONWRAPPER_HPP
#define PYTHONWRAPPER_HPP

// System includes
//
#include <Python.h>
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief This class is a small wrapper for the Python embedding
    */
   class PythonWrapper
   {
      public:
         /**
          * @brief Initialise the Python interpreter
          */
         static void init();

         /**
          * @brief Import the python module
          */
         static void import(const std::string& module);

         /**
          * @brief Set the python function object
          */
         static void setFunction(const std::string& func);

         /**
          * @brief Call the python function object with arguments
          */
         static PyObject* callFunction(PyObject* pArgs);

         /**
          * @brief Make a tuple from integer array
          */
         static PyObject* makeTuple(const ArrayI& val);

         /**
          * @brief Make a tuple from vector of double
          */
         static PyObject* makeTuple(const std::vector<MHDFloat>& val);

         /**
          * @brief Make a dictironary
          */
         static PyObject* makeDict(const std::vector<std::string>& key, const std::vector<MHDFloat>& val);

         /**
          * @brief Make a dictironary
          */
         static PyObject* makeDict(const std::map<std::string,int>& map);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> > &rList, PyObject *pList);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         static void fillMatrix(SparseMatrix& rMatrix, PyObject* pPyMat);

         /**
          * @brief Fill sparse matrix with data from Python call
          */
         static void fillMatrix(DecoupledZSparse& rMatrix, PyObject* pPyMat);
  
         /**
          * @brief Cleanup wrapper without finalize
          */
         static void cleanup();
  
         /**
          * @brief Finalise the Python interpreter
          */
         static void finalize();
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         PythonWrapper();

         /**
          * @brief Destructor
          */
         ~PythonWrapper();

         /**
          * @brief Python module object
          */
         static PyObject* mpModule;

         /**
          * @brief Python function object
          */
         static PyObject* mpFunc;
   };

}

#endif // PYTHONWRAPPER_HPP

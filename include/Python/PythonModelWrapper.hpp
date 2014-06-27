/**
 * @file PythonModelWrapper.hpp
 * @brief Static class wrapper to the Python model embedding
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PYTHONMODELWRAPPER_HPP
#define PYTHONMODELWRAPPER_HPP

// First include
//
#include "Python/PythonHeader.hpp"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Static class wrapper to the Python model embedding
    */
   class PythonModelWrapper
   {
      public:
         /**
          * @brief Initialise the Python interpreter
          */
         static void init();

         /**
          * @brief Import the Python module
          */
         static void import(const std::string& module);

         /**
          * @brief Create the Python model class
          */
         static void createModel(const std::string& model);

         /**
          * @brief Set the Python function object
          */
         static void setFunction(const std::string& func);

         /**
          * @brief Call the Python function object
          */
         static PyObject* callFunction();

         /**
          * @brief Call the Python function object with arguments
          */
         static PyObject* callFunction(PyObject* pArgs);

         /**
          * @brief Set the Python model method object
          */
         static void setMethod(const std::string& method);

         /**
          * @brief Call the Python model method object
          */
         static PyObject* callMethod();

         /**
          * @brief Call the Python model method object with arguments
          */
         static PyObject* callMethod(PyObject* pArgs);

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
         static void getList(std::vector<bool> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<NonDimensional::Id> &rList, PyObject *pList);

         /**
          * @brief Get data from list
          */
         static void getList(std::vector<PhysicalNames::Id> &rList, PyObject *pList);

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
         PythonModelWrapper();

         /**
          * @brief Destructor
          */
         ~PythonModelWrapper();

         /**
          * @brief Python module object
          */
         static PyObject* mpModule;

         /**
          * @brief Python function object
          */
         static PyObject* mpFunc;

         /**
          * @brief Python model object
          */
         static PyObject* mpModel;

         /**
          * @brief Python model method object
          */
         static PyObject* mpMethod;
   };

}

#endif // PYTHONMODELWRAPPER_HPP

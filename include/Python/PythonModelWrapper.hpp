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

namespace QuICC {

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
          * @brief Create the Python model class, first trying model+specialization
          */
         static void createModel(const std::string& model, const std::string& specialization);

         /**
          * @brief Set the Python function object
          */
         static void setFunction(const std::string& func);

         /**
          * @brief Set the Python function object
          */
         static void setFunction(const std::string& func, const std::string& submod);

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

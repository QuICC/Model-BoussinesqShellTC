/**
 * @file PythonWrapper.hpp
 * @brief Small wrapper for Python embedding
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PYTHONWRAPPER_HPP
#define PYTHONWRAPPER_HPP

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
          * @brief Import the Python module
          */
         static void import(const std::string& module);

         /**
          * @brief Create a Python class object
          */
         static void createClass(const std::string& name);

         /**
          * @brief Set the Python function object
          */
         static void setFunction(const std::string& func);

         /**
          * @brief Set the Python function from "sub module" object 
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
          * @brief Set the Python class method object
          */
         static void setMethod(const std::string& method);

         /**
          * @brief Call the Python class method object
          */
         static PyObject* callMethod();

         /**
          * @brief Call the Python class method object with arguments
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

         /**
          * @brief Python class object
          */
         static PyObject* mpClass;

         /**
          * @brief Python class method object
          */
         static PyObject* mpMethod;
   };

}

#endif // PYTHONWRAPPER_HPP

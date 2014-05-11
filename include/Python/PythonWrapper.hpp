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

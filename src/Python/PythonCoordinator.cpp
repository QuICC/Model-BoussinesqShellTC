/** 
 * @file PythonCoordinator.cpp
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
#include "Python/PythonCoordinator.hpp"

// Project includes
//
#include "Python/PythonConfig.hpp"
#include "Exceptions/Exception.hpp"
#include "IoTools/HumanToId.hpp"

namespace GeoMHDiSCC {

   int PythonCoordinator::sCounter = 0;

   void PythonCoordinator::init()
   {
      if(PythonCoordinator::sCounter == 0)
      {
         // Initialize the python interpreter
         Py_Initialize();

         // Setup the search path
         PyObject* sysPath = PySys_GetObject((char*)"path");
         PyList_Append(sysPath, PyUnicode_FromString(QUICC_PYTHON_DIR));
      }

      PythonCoordinator::registerWrapper();
   }

   void PythonCoordinator::registerWrapper()
   {
      PythonCoordinator::sCounter++;
   }

   void PythonCoordinator::unregisterWrapper()
   {
      PythonCoordinator::sCounter--;
   }

   void PythonCoordinator::finalize()
   {
      PythonCoordinator::unregisterWrapper();

      if(PythonCoordinator::sCounter == 0)
      {
         // Finalize
         Py_Finalize();
      }
   }

   PythonCoordinator::PythonCoordinator()
   {
   }

   PythonCoordinator::~PythonCoordinator()
   {
   }

}

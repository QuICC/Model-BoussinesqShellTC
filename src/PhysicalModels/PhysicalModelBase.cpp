/** 
 * @file PhysicalModelBase.cpp
 * @brief Source of the base for the physical models
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
#include "PhysicalModels/PhysicalModelBase.hpp"

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"
#include "IoTools/HumanToId.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> PhysicalModelBase::fieldIds(const std::string& pyName)
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(pyName);

      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      PythonWrapper::setFunction((char *)"all_fields");
      pValue = PythonWrapper::callFunction();

      // Create storage
      std::vector<PhysicalNames::Id> ids;
      PythonWrapper::getList(ids, pValue);

      // Clenup Python interpreter
      PythonWrapper::cleanup();

      return ids;
   }

   std::vector<NonDimensional::Id> PhysicalModelBase::paramIds(const std::string& pyName)
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(pyName);

      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      PythonWrapper::setFunction((char *)"nondimensional_parameters");
      pValue = PythonWrapper::callFunction();

      // Create storage
      std::vector<NonDimensional::Id> ids;
      PythonWrapper::getList(ids, pValue);

      // Cleanup Python interpreter
      PythonWrapper::cleanup();

      return ids;
   }

   std::vector<bool> PhysicalModelBase::isPeriodicBox(const std::string& pyName)
   {
      // Initialise Python interpreter
      PythonWrapper::init();

      // Load model module
      PythonWrapper::import(pyName);

      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      PythonWrapper::setFunction((char *)"periodicity");
      pValue = PythonWrapper::callFunction();

      // Create storage
      std::vector<bool> box;
      PythonWrapper::getList(box, pValue);

      // Cleanup Python interpreter
      PythonWrapper::cleanup();

      return box;
   }

}

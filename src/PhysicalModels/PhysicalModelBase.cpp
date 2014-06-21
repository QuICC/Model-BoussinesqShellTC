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
#include "Python/PythonModelWrapper.hpp"

namespace GeoMHDiSCC {

   std::vector<PhysicalNames::Id> PhysicalModelBase::fieldIds()
   {
      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      PythonModelWrapper::setMethod((char *)"all_fields");
      pValue = PythonModelWrapper::callMethod();

      // Create storage
      std::vector<PhysicalNames::Id> ids;
      PythonModelWrapper::getList(ids, pValue);

      // Clenup Python interpreter
      PythonModelWrapper::cleanup();

      return ids;
   }

   std::vector<NonDimensional::Id> PhysicalModelBase::paramIds()
   {
      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      PythonModelWrapper::setMethod((char *)"nondimensional_parameters");
      pValue = PythonModelWrapper::callMethod();

      // Create storage
      std::vector<NonDimensional::Id> ids;
      PythonModelWrapper::getList(ids, pValue);

      // Cleanup Python interpreter
      PythonModelWrapper::cleanup();

      return ids;
   }

   std::vector<bool> PhysicalModelBase::isPeriodicBox()
   {
      // Prepare Python call arguments
      PyObject *pValue;

      // Call model operator Python routine
      PythonModelWrapper::setMethod((char *)"periodicity");
      pValue = PythonModelWrapper::callMethod();

      // Create storage
      std::vector<bool> box;
      PythonModelWrapper::getList(box, pValue);

      // Cleanup Python interpreter
      PythonModelWrapper::cleanup();

      return box;
   }

}

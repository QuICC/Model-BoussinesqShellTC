/** \file CouplingInformation.cpp
 *  \brief Source of the base implementation of a scalar equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/CouplingInformation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   CouplingInformation::CouplingInformation()
      : mHasField(false)
   {
   }

   CouplingInformation::~CouplingInformation()
   {
   }

   void CouplingInformation::addField(const FieldComponents::Spectral::Id& component, const std::pair<PhysicalNames::Id, FieldComponents::Spectral::Id>& field)
   {
      // Set the coupling flag
      this->mHasField = true;

      // Store field coupling information
      this->mField.insert(std::make_pair(component,field));
   }

   void CouplingInformation::addInternal(const FieldComponents::Spectral::Id& component, const int nMat, const ArrayI& dim)
   {
      // Store internal coupling information
      this->mInternal.insert(std::make_pair(component,std::make_pair(nMat,dim)));
   }
}
}

/** 
 * @file IScalarEquation.cpp
 * @brief Source of the base implementation of a scalar equation
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
#include "Equations/IScalarEquation.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   IScalarEquation::IScalarEquation(SharedEquationParameters spEqParams)
      : IEquation(spEqParams)
   {
   }

   IScalarEquation::~IScalarEquation()
   {
   }

   void IScalarEquation::setUnknown(Datatypes::SharedScalarVariableType spUnknown)
   {
      this->mspUnknown = spUnknown;
   }

   const Datatypes::ScalarVariableType& IScalarEquation::unknown() const
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   Datatypes::ScalarVariableType& IScalarEquation::rUnknown()
   {
      // Safety assert
      assert(this->mspUnknown);

      return *this->mspUnknown;
   }

   SharedResolution IScalarEquation::spRes() const
   {
      return this->unknown().dom(0).spRes();
   }

   void IScalarEquation::updateDealiasedUnknown(const Datatypes::SpectralScalarType& rhs, FieldComponents::Spectral::Id compId)
   {
      // Assert scalar
      assert(compId == FieldComponents::Spectral::SCALAR);
      // Assert dealiasing has taken place!
      assert(this->rUnknown().rDom(0).rPerturbation().data().rows() < rhs.data().rows());
      assert(this->rUnknown().rDom(0).rPerturbation().data().cols() == rhs.data().cols());

      // Copy values over into unknown
      this->rUnknown().rDom(0).rPerturbation().setData(rhs.data().topRows(this->rUnknown().rDom(0).rPerturbation().data().rows()));
   }

   void IScalarEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      // Store the boundary condition list
      this->mspBcIds = spBcIds;

      // Make sure it is safe to do nothing
      bool needInit = this->couplingInfo(FieldComponents::Spectral::SCALAR).hasQuasiInverse();

      CouplingInformation::FieldId_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).explicitRange();
      needInit = needInit && (fRange.first == fRange.second);

      // Initialise spectral matrices
      if(needInit)
      {
         if(EquationToolsType::EIGEN_DIMS == 0)
         {
            this->initSpectralMatricesNoEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else if(EquationToolsType::EIGEN_DIMS == 1)
         {
            this->initSpectralMatrices1DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else if(EquationToolsType::EIGEN_DIMS == 2)
         {
            this->initSpectralMatrices2DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else if(EquationToolsType::EIGEN_DIMS == 3)
         {
            this->initSpectralMatrices3DEigen(spBcIds, FieldComponents::Spectral::SCALAR);
         } else
         {
            throw Exception("initSpectralMatrices: unhandled number of eigen dimensions!");
         }
      }
   }
}
}

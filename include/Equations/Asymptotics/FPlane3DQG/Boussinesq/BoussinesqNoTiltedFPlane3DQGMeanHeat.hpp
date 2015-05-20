/**
 * @file BoussinesqNoTiltedFPlane3DQGMeanHeat.hpp
 * @brief Implementation of the mean heat computation for the Boussinesq tilted F-plane 3DQG model 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQNOTILTEDFPLANE3DQGMEANHEAT_HPP
#define BOUSSINESQNOTILTEDFPLANE3DQGMEANHEAT_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the mean heat computation for the Boussinesq tilted F-plane 3DQG model
    */
   class BoussinesqNoTiltedFPlane3DQGMeanHeat: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         BoussinesqNoTiltedFPlane3DQGMeanHeat(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqNoTiltedFPlane3DQGMeanHeat();
         
         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param id      ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const;

         /**
          * @brief Compute the source term
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         virtual Datatypes::SpectralScalarType::PointType sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

      private:
   };

}
}

#endif // BOUSSINESQNOTILTEDFPLANE3DQGMEANHEAT_HPP

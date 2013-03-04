/** \file Beta3DQGStreamfunction.hpp
 *  \brief Implementation of the streamfunction equation for the 3DQG beta model
 */

#ifndef BETA3DQGSTREAMFUNCTION_HPP
#define BETA3DQGSTREAMFUNCTION_HPP

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
#include "Equations/IScalarEquation.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * \brief Implementation of the streamfunction equation for the 3DQG beta model
    */
   class Beta3DQGStreamfunction: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Beta3DQGStreamfunction(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Beta3DQGStreamfunction();

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const;

         /**
          * @brief Compute the linear term
          *
          * @param rRHS    RHS of timestepping equation
          */
         virtual void computeLinear(Datatypes::SpectralScalarType& rRHS) const;

         /**
          * @brief Transfer timestepper output to unknown and update vorticity
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Transfer timestepper output to unknown and update vorticity
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Set equation matrices
          */
         virtual void setSpectralMatrices(Spectral::SpectralSelector<Dimensions::Transform::TRA1D>::Type& spec1D, Spectral::SpectralSelector<Dimensions::Transform::TRA2D>::Type& spec2D, Spectral::SpectralSelector<Dimensions::Transform::TRA3D>::Type& spec3D);
         
      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
   };

}
}

#endif // BETA3DQGSTREAMFUNCTION_HPP

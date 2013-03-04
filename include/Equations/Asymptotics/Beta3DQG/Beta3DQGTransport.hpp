/** \file Beta3DQGTransport.hpp
 *  \brief Implementation of the transport equation for the 3DQG beta model
 */

#ifndef BETA3DQGTRANSPORT_HPP
#define BETA3DQGTRANSPORT_HPP

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
    * \brief Implementation of the transport equation for the 3DQG beta model
    */
   class Beta3DQGTransport: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Beta3DQGTransport(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Beta3DQGTransport();

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
          * @brief Prepare the RHS for the timestep computation
          *
          * @param rhs    RHS of timestepping equation
          */
         virtual void prepareTimestep(const Datatypes::SpectralScalarType& rhs);

         /**
          * @brief Method to transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Method to transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Set the equation matrices
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

#endif // BETA3DQGTRANSPORT_HPP

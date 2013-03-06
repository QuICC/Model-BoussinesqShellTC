/** \file IEvolutionEquation.hpp
 *  \brief Base building block for the implementation of a time dependend evolution equation
 *
 *  \mhdBug Needs test
 */

#ifndef IEVOLUTIONEQUATION_HPP
#define IEVOLUTIONEQUATION_HPP

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
#include "Enums/Dimensions.hpp"
#include "Enums/PhysicalNames.hpp"
#include "Enums/FieldComponents.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"
#include "Equations/IEquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "Equations/EquationData.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Variables/VariableRequirement.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of a time dependend evolution equation
    */
   class IEvolutionEquation : public EquationData
   {
      public:
         /// Typedef for the boundary conditions map key type
         typedef std::pair<FieldComponents::Spectral::Id, Dimensions::Simulation::Id> BcKeyType;

         /// Typedef for the boundary conditions map storage type
         typedef std::vector<std::pair<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position> > BcMapType;

         /// Typedef for the boundary conditions map storage type
         typedef std::map<BcKeyType, BcMapType> BcEqMapType;

         /**
          * @brief Simple constructor
          */
         explicit IEvolutionEquation(SharedIEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~IEvolutionEquation();

         /**
          * @brief Initialise the equation
          */
         virtual void init();

         /**
          * @brief Pure virtual method to transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Pure virtual method to transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Pure virtual method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Pure virtual method to transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Set the equation matrices
          *
          * @param bcIds   List of boundary condition IDs
          * @param cbcIds  List of coupled boundary condition IDs
          */
         virtual void setSpectralMatrices(const BcEqMapType& bcIds, const std::map<PhysicalNames::Id, BcEqMapType>& cbcIds) = 0;

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Transfer equation input to timestepper
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTInput(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Transfer timestepper output to equation unknown
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void copyTOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start) = 0;
         
      protected:
         /**
          * @brief Apply quasi-inverse to nonlinear values
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void applyQuasiInverse(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Apply quasi-inverse to nonlinear values
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start indx for the storage
          */
         virtual void applyQuasiInverse(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Set the equation variable requirements
          */
         virtual void setRequirements() = 0;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() = 0;

      private:
   };

   /// Typedef for a smart IEvolutionEquation
   typedef SharedPtrMacro<IEvolutionEquation> SharedIEvolutionEquation;
}
}

#endif // IEVOLUTIONEQUATION_HPP

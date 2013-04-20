/** \file IEvolutionEquation.hpp
 *  \brief Base building block for the implementation of a time dependend evolution equation
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
#include "Enums/FieldIds.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"
#include "Equations/IEquationParameters.hpp"
#include "Equations/CouplingInformation.hpp"
#include "Equations/EquationData.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Variables/VariableRequirement.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Base building block for the implementation of a time dependend evolution equation
    */
   class IEvolutionEquation : public EquationData
   {
      public:
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
          * @brief Compute the explicit linear terms
          *
          * @param id            Component ID
          * @param storage       Storage for the equation values
          * @param linearField   Field with explicit linear dependence 
          * @param matIdx        Index of the given data
          * @param start         Start index for the storage
          */
         virtual void computeLinear(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const DecoupledZMatrix& linearField, const int matIdx, const int start) const;

         /**
          * @brief Compute the explicit linear terms
          *
          * @param id            Component ID
          * @param storage       Storage for the equation values
          * @param linearField   Field with explicit linear dependence 
          * @param matIdx        Index of the given data
          * @param start         Start index for the storage
          */
         virtual void computeLinear(FieldComponents::Spectral::Id id, MatrixZ& storage, const MatrixZ& linearField, const int matIdx, const int start) const;

         /**
          * @brief Pure virtual method to transfer equation input to timestepper
          *
          * @param id         Component ID
          * @param storage    Storage for the equation values
          * @param matIdx     Index of the given data
          * @param start      Start indx for the storage
          */
         virtual void timestepInput(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start) = 0;

         /**
          * @brief Pure virtual method to transfer equation input to timestepper
          *
          * @param id         Component ID
          * @param storage    Storage for the equation values
          * @param matIdx     Index of the given data
          * @param start      Start indx for the storage
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

         /**
          * @brief Apply quasi-inverse to the nonlinear term
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void applyNLQuasiInverse(FieldComponents::Spectral::Id id, DecoupledZMatrix& storage, const int matIdx, const int start);

         /**
          * @brief Apply quasi-inverse to the nonlinear term
          *
          * @param id      Component ID
          * @param storage Storage for the equation values
          * @param matIdx  Index of the given data
          * @param start   Start index for the storage
          */
         virtual void applyNLQuasiInverse(FieldComponents::Spectral::Id id, MatrixZ& storage, const int matIdx, const int start);

         /**
          * @brief Get linear operator row of coupled matrix
          */
         virtual DecoupledZSparse linearRow(FieldComponents::Spectral::Id id, const int matIdx) const = 0;

         /**
          * @brief Get time operator row of coupled matrix
          */
         virtual DecoupledZSparse timeRow(FieldComponents::Spectral::Id id, const int matIdx) const = 0;

         /**
          * @brief Get boundary operator row of coupled matrix
          */
         virtual DecoupledZSparse boundaryRow(FieldComponents::Spectral::Id id, const int matIdx) const = 0;

         /**
          * @brief Initialise the spectral equation matrices
          *
          * @param spBcIds   List of boundary condition IDs
          */
         virtual void initSpectralMatrices(const SharedSimulationBoundary spBcIds) = 0;
         
      protected:
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

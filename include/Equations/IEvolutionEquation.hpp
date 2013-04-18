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
#include "Enums/PhysicalNames.hpp"
#include "Enums/FieldComponents.hpp"
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
          */
         virtual void setSpectralMatrices(const SimulationBoundary& bcIds) = 0;

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

         /**
          * @brief Is coupled system of equation complex?
          */
         virtual bool isSystemComplex(FieldComponents::Spectral::Id id) const = 0;

         /**
          * @brief Get starting index for equation (useful to exclude m=0 for example)
          */
         virtual int startIndex(FieldComponents::Spectral::Id id) const = 0;

         /**
          * @brief Get index of the system of equations
          */
         virtual int systemIndex(FieldComponents::Spectral::Id id) const = 0;

         /**
          * @brief Get size of the system of equations
          */
         virtual int systemSize(FieldComponents::Spectral::Id id) const = 0;

         /**
          * @brief Get linear operator row of coupled matrix
          */
         virtual linearRow(FieldComponents::Spectral::Id id, const int matIdx) const = 0;

         /**
          * @brief Get time operator row of coupled matrix
          */
         virtual timeRow(FieldComponents::Spectral::Id id, const int matIdx) const = 0;

         /**
          * @brief Get boundary operator row of coupled matrix
          */
         virtual boundaryRow(FieldComponents::Spectral::Id id, const int matIdx) const = 0;


      private:
   };

   /// Typedef for a smart IEvolutionEquation
   typedef SharedPtrMacro<IEvolutionEquation> SharedIEvolutionEquation;
}
}

#endif // IEVOLUTIONEQUATION_HPP

/** \file SimulationBoundary.hpp
 *  \brief Implementation of a simple simulation wide boundary condition interface
 */

#ifndef SIMULATIONBOUNDARY_HPP
#define SIMULATIONBOUNDARY_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Runtime.hpp"
#include "IoControl/ControlInterface.hpp"
#include "Equations/IEvolutionEquation.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a simple simulation wide boundary condition interface
    */
   class SimulationBoundary
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationBoundary();

         /**
          * @brief Destructor
          */
         ~SimulationBoundary();

         /**
          * @brief Get boundary conditions for specific equation
          */
         const Equations::IEvolutionEquation::BcEqMapType& bc(const PhysicalNames::Id id) const;

         /**
          * @brief Get coupled boundary conditions for specific equation
          */
         const std::map<PhysicalNames::Id, Equations::IEvolutionEquation::BcEqMapType>& cbc(const PhysicalNames::Id id) const;

         /**
          * @brief Initialise storage for a new equation
          *
          * @param id ID of the equation variable
          */
         void initStorage(const PhysicalNames::Id id);

         /**
          * @brief Initialise storage for a new boundary condition component
          *
          * @param id ID of the equation variable
          */
         void initBcStorage(const PhysicalNames::Id id, const Equations::IEvolutionEquation::BcKeyType& key);

         /**
          * @brief Initialise storage for a new couple boundary condition component
          *
          * @param eqId ID of the equation variable
          * @param cpId ID of the coupled variable
          */
         void initCBcStorage(const PhysicalNames::Id eqId, const PhysicalNames::Id cpId, const Equations::IEvolutionEquation::BcKeyType& key);

         /**
          * @brief Add boundary condition
          *
          * @param id ID of the equation variable
          */
         void addBc(const PhysicalNames::Id id, const Equations::IEvolutionEquation::BcKeyType& key, Spectral::BoundaryConditions::Id bcId, Spectral::IBoundary::Position pos);

         /**
          * @brief Add coupled condition
          *
          * @param eqId ID of the equation variable
          * @param cpId ID of the coupled variable
          */
         void addCBc(const PhysicalNames::Id eqId, const PhysicalNames::Id cpId, const Equations::IEvolutionEquation::BcKeyType& key, Spectral::BoundaryConditions::Id bcId, Spectral::IBoundary::Position pos);
         
      protected:

      private:
         /**
          * @brief Storage for the boundary conditions
          */
         std::map<PhysicalNames::Id, Equations::IEvolutionEquation::BcEqMapType> mBcs;

         /**
          * @brief Storage for the coupled boundary conditions
          */
         std::map<PhysicalNames::Id, std::map<PhysicalNames::Id, Equations::IEvolutionEquation::BcEqMapType> > mCBcs;
   };
}

#endif // SIMULATIONBOUNDARY_HPP

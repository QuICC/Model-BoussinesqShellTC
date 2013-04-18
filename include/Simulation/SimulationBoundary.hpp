/** \file SimulationBoundary.hpp
 *  \brief Implementation of a simple simulation wide boundary condition interface
 */

#ifndef SIMULATIONBOUNDARY_HPP
#define SIMULATIONBOUNDARY_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldComponents.hpp"
#include "Enums/PhysicalNames.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a simple simulation wide boundary condition interface
    */
   class SimulationBoundary
   {
      public:
         /// Typedef for the boundary conditions map key type
         typedef std::pair<FieldComponents::Spectral::Id, Dimensions::Simulation::Id> BcKeyType;

         /// Typedef for the boundary conditions map storage type
         typedef std::vector<std::pair<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position> > BcMapType;

         /// Typedef for the boundary conditions map storage type
         typedef std::map<BcKeyType, BcMapType> BcEqMapType;

         /**
          * @brief Constructor
          */
         SimulationBoundary();

         /**
          * @brief Destructor
          */
         ~SimulationBoundary();

         /**
          * @brief Check if boundary condition for equation exist
          *
          * @param eqId    Physical ID for equation
          */
         bool hasEquation(const PhysicalNames::Id eqId) const;

         /**
          * @brief Check if boundary condition for field in given equation exist
          *
          * @param eqId       Physical ID for equation
          * @param fieldId    Physical ID for field
          */
         bool hasField(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId) const;

         /**
          * @brief Get boundary conditions for specific equation
          *
          * @param eqId       Physical ID for equation
          * @param fieldId    Physical ID for field
          */
         const SimulationBoundary::BcEqMapType& bcs(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId) const;

         /**
          * @brief Initialise storage for a new equation
          *
          * @param id ID of the equation
          */
         void initStorage(const PhysicalNames::Id id);

         /**
          * @brief Initialise storage for a new boundary condition component
          *
          * @param eqId    ID of the equation
          * @param fieldId ID of the variable
          */
         void initBcStorage(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId, const SimulationBoundary::BcKeyType& key);

         /**
          * @brief Add boundary condition
          *
          * @param eqId    ID of the equation
          * @param fieldId ID of the variable
          */
         void addBc(const PhysicalNames::Id eqId, const PhysicalNames::Id fieldId, const SimulationBoundary::BcKeyType& key, Spectral::BoundaryConditions::Id bcId, Spectral::IBoundary::Position pos);
         
      protected:

      private:
         /**
          * @brief Storage for the coupled boundary conditions
          */
         std::map<PhysicalNames::Id, std::map<PhysicalNames::Id, SimulationBoundary::BcEqMapType> > mBcs;
   };
}

#endif // SIMULATIONBOUNDARY_HPP

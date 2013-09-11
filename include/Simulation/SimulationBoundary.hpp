/** 
 * @file SimulationBoundary.hpp
 * @brief Implementation of a simple simulation wide boundary condition interface
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "Enums/FieldIds.hpp"
#include "SpectralOperators/BoundaryConditions.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a simple simulation wide boundary condition interface
    */
   class SimulationBoundary
   {
      public:
         /// Typedef for the boundary conditions map storage type
         typedef std::vector<std::pair<Spectral::BoundaryConditions::Id,Spectral::IBoundary::Position> > BcMapType;

         /// Typedef for the boundary conditions map storage type
         typedef std::map<Dimensions::Simulation::Id, BcMapType> BcEqMapType;

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
         bool hasEquation(const SpectralFieldId eqId) const;

         /**
          * @brief Check if boundary condition for field in given equation exist
          *
          * @param eqId       Physical ID for equation
          * @param fieldId    Physical ID for field
          */
         bool hasField(const SpectralFieldId eqId, const SpectralFieldId fieldId) const;

         /**
          * @brief Get boundary conditions for specific equation
          *
          * @param eqId       Physical ID for equation
          * @param fieldId    Physical ID for field
          */
         const SimulationBoundary::BcEqMapType& bcs(const SpectralFieldId eqId, const SpectralFieldId fieldId) const;

         /**
          * @brief Initialise storage for a new equation
          *
          * @param id ID of the equation
          */
         void initStorage(const SpectralFieldId id);

         /**
          * @brief Initialise storage for a new boundary condition component
          *
          * @param eqId    ID of the equation
          * @param fieldId ID of the variable
          * @param dimId   ID of the dimension
          */
         void initBcStorage(const SpectralFieldId eqId, const SpectralFieldId fieldId, const Dimensions::Simulation::Id dimId);

         /**
          * @brief Add boundary condition
          *
          * @param eqId    ID of the equation
          * @param fieldId ID of the variable
          * @param dimId   ID of the dimension
          */
         void addBc(const SpectralFieldId eqId, const SpectralFieldId fieldId, const Dimensions::Simulation::Id dimId, Spectral::BoundaryConditions::Id bcId, Spectral::IBoundary::Position pos);
         
      protected:

      private:
         /**
          * @brief Storage for the coupled boundary conditions
          */
         std::map<SpectralFieldId, std::map<SpectralFieldId, SimulationBoundary::BcEqMapType> > mBcs;
   };

   /// Typedef for a shared pointer to a SimulationBoundary object
   typedef SharedPtrMacro<SimulationBoundary>   SharedSimulationBoundary;
}

#endif // SIMULATIONBOUNDARY_HPP

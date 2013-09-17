/**
 * @file BoussinesqPerBetaCylGSystem.hpp
 * @brief Implementation of the matrix blocks for the Boussinesq beta model with cylindrical gravity with periodic radius 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BOUSSINESQPERBETACYLGSYSTEM_HPP
#define BOUSSINESQPERBETACYLGSYSTEM_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Equations/CouplingInformation.hpp"
#include "Simulation/SimulationBoundary.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the matrix blocks for the Boussinesq beta model with cylindrical gravity with periodic radius
    *
    * The implementation of all operators assume that all equations are written in the form
    * \f$ L - T = NL \f$
    * where L are the linear operators, T the time derivatives and NL the nonlinear terms
    */
   class BoussinesqPerBetaCylGSystem
   {
      public:
         /**
          * @brief Set the coupling information
          *
          * @param rInfo   Output coupling information
          * @param eqId    Physical ID of the equation
          * @param nZ      Matrix size in X
          * @param nX      Number of radial wave numbers
          * @param nY      Number of azimuthal wave numbers
          */
         static void setCouplingInfo(CouplingInformation& rInfo, const SpectralFieldId eqId, const int nZ, const int modes);

         /**
          * @brief Get the quasi-inverse matrix operator for an equation
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param nZ      Matrix size in Z
          */
         static void quasiInverse(SparseMatrix& mat, const SpectralFieldId eqId, const int nZ);

         /**
          * @brief Get the time matrix block for an equation
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param nZ      Matrix size in Z
          * @param kX      Radial wave number
          * @param kY      Azimuthal wave number
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void timeBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const int nZ, const MHDFloat kX, const MHDFloat kY, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);

         /**
          * @brief Get the linear matrix block for an equation on given field
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param fieldId Physical ID of the field
          * @param nZ      Matrix size in Z
          * @param kX      Radial wave number
          * @param kY      Azimuthal wave number
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void linearBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const int nZ, const MHDFloat kX, const MHDFloat kY, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);

         /**
          * @brief Get the boundary condition matrix block for an equation on given field
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param fieldId Physical ID of the field
          * @param spBcIds Boundary condition IDs
          * @param nZ      Matrix size in Z
          * @param kX      Radial wave number
          * @param kY      Azimuthal wave number
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void boundaryBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nZ, const MHDFloat kX, const MHDFloat kZ, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);
         
      protected:

      private:
         /**
          * @brief Simple constructor
          */
         BoussinesqPerBetaCylGSystem();

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqPerBetaCylGSystem();
   };

}
}

#endif // BOUSSINESQPERBETACYLGSYSTEM_HPP

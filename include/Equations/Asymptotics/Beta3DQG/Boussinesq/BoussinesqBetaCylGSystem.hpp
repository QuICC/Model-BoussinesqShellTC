/** \file BoussinesqBetaCylGSystem.hpp
 *  \brief Implementation of the matrix blocks for the 3DQG beta model
 */

#ifndef BETA3DQGSYSTEM_HPP
#define BETA3DQGSYSTEM_HPP

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
    * \brief Implementation of the matrix blocks for the 3DQG beta model
    *
    * The implementation of all operators assume that all equations are written in the form
    * \f$ L - T = NL \f$
    * where L are the linear operators, T the time derivatives and NL the nonlinear terms
    */
   class BoussinesqBetaCylGSystem
   {
      public:
         /**
          * @brief Set the coupling information
          *
          * @param rInfo   Output coupling information
          * @param eqId    Physical ID of the equation
          * @param nx      Matrix size in X
          * @param nz      Matrix size in Z
          * @param ny      Number of azimuthal wave numbers
          */
         static void setCouplingInfo(CouplingInformation& rInfo, const SpectralFieldId eqId, const int nx, const int nz, const int ny);

         /**
          * @brief Get the quasi-inverse matrix operator for an equation
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param nx      Matrix size in X
          * @param nz      Matrix size in Z
          */
         static void quasiInverse(SparseMatrix& mat, const SpectralFieldId eqId, const int nx, const int nz);

         /**
          * @brief Get the time matrix block for an equation
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param nx      Matrix size in X
          * @param nz      Matrix size in Z
          * @param k       Wave number k
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void timeBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const int nx, const int nz, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);

         /**
          * @brief Get the linear matrix block for an equation on given field
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param fieldId Physical ID of the field
          * @param nx      Matrix size in X
          * @param nz      Matrix size in Z
          * @param k       Wave number k
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void linearBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const int nx, const int nz, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);

         /**
          * @brief Get the boundary condition matrix block for an equation on given field
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param fieldId Physical ID of the field
          * @param spBcIds Boundary condition IDs
          * @param nx      Matrix size in X
          * @param nz      Matrix size in Z
          * @param k       Wave number k
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void boundaryBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nx, const int nz, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);
         
      protected:

      private:
         /**
          * @brief Simple constructor
          */
         BoussinesqBetaCylGSystem();

         /**
          * @brief Simple empty destructor
          */
         virtual ~BoussinesqBetaCylGSystem();
   };

}
}

#endif // BETA3DQGSYSTEM_HPP

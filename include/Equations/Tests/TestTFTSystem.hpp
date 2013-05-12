/** \file TestTFTSystem.hpp
 *  \brief Implementation of the matrix blocks for the test equations for the TFT scheme
 */

#ifndef TESTTFTSYSTEM_HPP
#define TESTTFTSYSTEM_HPP

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
    * @brief Implementation of the matrix blocks for the test equations for the TFT scheme
    *
    * The implementation of all operators assume that all equations are written in the form
    * \f$ L - T = NL \f$
    * where L are the linear operators, T the time derivatives and NL the nonlinear terms
    */
   class TestTFTSystem
   {
      public:
         /**
          * @brief Set the coupling information
          *
          * @param rInfo   Output coupling information
          * @param eqId    Physical ID of the equation
          * @param nX      Matrix size in X
          * @param nZ      Matrix size in Z
          * @param nY      Number of azimuthal wave numbers
          */
         static void setCouplingInfo(CouplingInformation& rInfo, const SpectralFieldId eqId, const int nX, const int nZ, const int nY);

         /**
          * @brief Get the quasi-inverse matrix operator for an equation
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param nX      Matrix size in X
          * @param nZ      Matrix size in Z
          */
         static void quasiInverse(SparseMatrix& mat, const SpectralFieldId eqId, const int nX, const int nZ);

         /**
          * @brief Get the time matrix block for an equation
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param nX      Matrix size in X
          * @param nZ      Matrix size in Z
          * @param k       Wave number k
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void timeBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const int nX, const int nZ, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);

         /**
          * @brief Get the linear matrix block for an equation on given field
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param fieldId Physical ID of the field
          * @param nX      Matrix size in X
          * @param nZ      Matrix size in Z
          * @param k       Wave number k
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void linearBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const int nX, const int nZ, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);

         /**
          * @brief Get the boundary condition matrix block for an equation on given field
          *
          * @param mat     Storage for output matrix
          * @param eqId    Physical ID of the equation
          * @param fieldId Physical ID of the field
          * @param spBcIds Boundary condition IDs
          * @param nX      Matrix size in X
          * @param nZ      Matrix size in Z
          * @param k       Wave number k
          * @param Ra      Rayleigh number
          * @param Pr      Prandtl number
          * @param Gamma   Gamma number
          * @param chi     Angle Chi
          */
         static void boundaryBlock(DecoupledZSparse& mat, const SpectralFieldId eqId, const SpectralFieldId fieldId, const SharedSimulationBoundary spBcIds, const int nX, const int nZ, const MHDFloat k, const MHDFloat Ra, const MHDFloat Pr, const MHDFloat Gamma, const MHDFloat chi);
         
      protected:

      private:
         /**
          * @brief Simple constructor
          */
         TestTFTSystem();

         /**
          * @brief Simple empty destructor
          */
         virtual ~TestTFTSystem();
   };

}
}

#endif // TESTTFTSYSTEM_HPP

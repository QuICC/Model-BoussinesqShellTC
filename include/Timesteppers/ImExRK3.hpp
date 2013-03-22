/** \file ImExRK3.hpp
 *  \brief Implementation of an implicit/explicit Runge-Kutta scheme of order ~3
 */

#ifndef IMEXRK3_HPP
#define IMEXRK3_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   /**
    * \brief Implementation of an implicit/explicit Runge-Kutta scheme of order ~3
    */
   class ImExRK3
   {
      public:
         /**
          * @brief Get factor for mass matrix on LHS
          *
          * @param step Substep
          */
         static MHDFloat lhsT(const int step);

         /**
          * @brief Get factor for linear matrix on LHS
          *
          * @param step Substep
          */
         static MHDFloat lhsL(const int step);

         /**
          * @brief Get factor for mass matrix on RHS
          *
          * @param step Substep
          */
         static MHDFloat rhsT(const int step);

         /**
          * @brief Get factor for linear matrix on RHS
          *
          * @param step Substep
          */
         static MHDFloat rhsL(const int step);

         /**
          * @brief Get factor for nonlinear term
          *
          * @param step Substep
          */
         static MHDFloat rhsN(const int step);

         /**
          * @brief Get factor for previous nonlinear term
          *
          * @param step Substep
          */
         static MHDFloat rhsNN(const int step);

         /**
          * @brief Number of substeps for final step
          */
         static const int STEPS;
         
      protected:
         /**
          * @brief Storage for the alpha parameters
          */
         static const Eigen::Array<MHDFloat,3,1> mAlpha;

         /**
          * @brief Storage for the beta parameters
          */
         static const Eigen::Array<MHDFloat,3,1> mBeta;

         /**
          * @brief Storage for the gamma parameters
          */
         static const Eigen::Array<MHDFloat,3,1> mGamma;

         /**
          * @brief Storage for the zeta parameters
          */
         static const Eigen::Array<MHDFloat,3,1> mZeta;

      private:
         /**
          * @brief Constructor
          */
         ImExRK3();

         /**
          * @brief Destructor
          */
         ~ImExRK3();

   };
}
}

#endif // IMEXRK3_HPP

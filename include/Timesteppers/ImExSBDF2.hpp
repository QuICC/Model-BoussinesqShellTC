/** 
 * @file ImExSBDF2.hpp
 * @brief Implementation of an implicit/explicit SBDF scheme of order 2
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IMEXSBDF2_HPP
#define IMEXSBDF2_HPP

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
    * @brief Implementation of an implicit/explicit SBDF scheme of order 2
    */
   class ImExSBDF2
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
          * @brief Get factor for mass matrix on RHS at t_(n-i)
          *
          * @param time Time index
          * @param step Substep
          */
         static MHDFloat rhsT(const int i, const int step);

         /**
          * @brief Get factor for linear matrix on RHS at t_(n-i)
          *
          * @param time Time index
          * @param step Substep
          */
         static MHDFloat rhsL(const int i, const int step);

         /**
          * @brief Get factor for nonlinear term at t_(n-i)
          *
          * @param i    Time index
          * @param step Substep
          */
         static MHDFloat rhsN(const int i, const int step);

         /**
          * @brief Field memory needed at current step
          *
          * @param step Substep
          */
         static int fieldMemory(const int step);

         /**
          * @brief Nonlinear term memory needed at current step
          *
          * @param step Substep
          */
         static int nonlinearMemory(const int step);

         /**
          * @brief Number of substeps for final step
          */
         static const int STEPS;

         /**
          * @brief Number of previous field values required
          */
         static const int FIELD_MEMORY;

         /**
          * @brief Number of previous nonlinear terms required
          */
         static const int NONLINEAR_MEMORY;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         ImExSBDF2();

         /**
          * @brief Destructor
          */
         ~ImExSBDF2();

   };
}
}

#endif // IMEXSBDF2_HPP

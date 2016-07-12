/** 
 * @file ImExSBDF2.cpp
 * @brief Implementation of an implicit/explicit SBDF scheme of order 2
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExSBDF2.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Timestep {

   // Scheme requires single step
   const int ImExSBDF2::STEPS = 1;

   // Scheme order
   const int ImExSBDF2::ORDER = 2;

   // Scheme requires field at t_(n-1)
   const int ImExSBDF2::FIELD_MEMORY = 1;

   // Scheme requires nonlinear term at t(n-1)
   const int ImExSBDF2::NONLINEAR_MEMORY = 1;

   // Name of the scheme
   const std::string ImExSBDF2::NAME = "ImExSBDF2";

   MHDFloat ImExSBDF2::lhsT(const int step)
   {
      assert(step < ImExSBDF2::STEPS);

      return 3.0/2.0;
   }

   MHDFloat ImExSBDF2::lhsL(const int step)
   {
      assert(step < ImExSBDF2::STEPS);

      return 1.0;
   }

   MHDFloat ImExSBDF2::rhsT(const int i, const int step)
   {
      assert(step < ImExSBDF2::STEPS);
      assert(i > -1);
      assert(i < ImExSBDF2::FIELD_MEMORY+1);

      MHDFloat coeff;
      if(i == 0)
      {
         coeff = 2.0;
      } else
      {
         coeff = -1.0/2.0;
      }

      return coeff;
   }

   MHDFloat ImExSBDF2::rhsL(const int i, const int step)
   {
      assert(step < ImExSBDF2::STEPS);
      assert(i > -1);
      assert(i < ImExSBDF2::FIELD_MEMORY+1);

      return 0.0;
   }

   MHDFloat ImExSBDF2::rhsN(const int i, const int step)
   {
      assert(step < ImExSBDF2::STEPS);
      assert(i > -1);
      assert(i < ImExSBDF2::NONLINEAR_MEMORY+1);

      MHDFloat coeff;
      if(i == 0)
      {
         coeff = 2.0;
      } else
      {
         coeff = -1.0;
      }

      return coeff;
   }

   int ImExSBDF2::fieldMemory(const int step)
   {
      assert(step < ImExSBDF2::STEPS);

      return ImExSBDF2::FIELD_MEMORY;
   }

   int ImExSBDF2::nonlinearMemory(const int step)
   {
      assert(step < ImExSBDF2::STEPS);

      return ImExSBDF2::NONLINEAR_MEMORY;
   }

   void ImExSBDF2::init()
   {
   }

   void ImExSBDF2::useEmbedded()
   {
      throw Exception("Tried to activate inexistant embedded scheme!");
   }

}
}

/** 
 * @file ImExRK2CN.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta 2nd order and Crank-Nicolson scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExRK2CN.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Timestep {

   // Scheme requires 3 substeps
   const int ImExRK2CN::STEPS = 2;

   // Scheme order
   const int ImExRK2CN::ORDER = 2;

   // Scheme has embedded lower order scheme?
   const bool ImExRK2CN::HAS_EMBEDDED = false;

   // Name of the scheme
   const std::string ImExRK2CN::NAME = "ImExRK2CN";

   // Set the implicit a factors
   Eigen::Array<MHDFloat,2,2> ImExRK2CN::mAIm = Eigen::Array<MHDFloat,2,2>::Zero();

   // Set the implicit b factors
   Eigen::Array<MHDFloat,2,1> ImExRK2CN::mBIm = Eigen::Array<MHDFloat,2,1>::Zero();

   // Set the explicit a factors
   Eigen::Array<MHDFloat,2,2> ImExRK2CN::mAEx = Eigen::Array<MHDFloat,2,2>::Zero();

   // Set the explicit b factors
   Eigen::Array<MHDFloat,2,1> ImExRK2CN::mBEx = Eigen::Array<MHDFloat,2,1>::Zero();

   // Set the implicit embedded b factors
   Eigen::Array<MHDFloat,2,1> ImExRK2CN::mBImErr = Eigen::Array<MHDFloat,2,1>::Zero();

   // Set the explicit embedded b factors
   Eigen::Array<MHDFloat,2,1> ImExRK2CN::mBExErr = Eigen::Array<MHDFloat,2,1>::Zero();

   MHDFloat ImExRK2CN::aIm(const int i, const int j)
   {
      return ImExRK2CN::mAIm(i,j);
   }

   MHDFloat ImExRK2CN::bIm(const int i)
   {
      return ImExRK2CN::mBIm(i);
   }

   MHDFloat ImExRK2CN::aEx(const int i, const int j)
   {
      return ImExRK2CN::mAEx(i,j);
   }

   MHDFloat ImExRK2CN::bEx(const int i)
   {
      return ImExRK2CN::mBEx(i);
   }

   MHDFloat ImExRK2CN::bImErr(const int i)
   {
      return ImExRK2CN::mBImErr(i);
   }

   MHDFloat ImExRK2CN::bExErr(const int i)
   {
      return ImExRK2CN::mBExErr(i);
   }

   void ImExRK2CN::init()
   {
      // Initialize implicit a factors
      ImExRK2CN::mAIm(1,1) = 2./5.;

      // Initialize implicit b factors
      ImExRK2CN::mBIm(1) = 5./6.;

      // Initialize explicit a factors
      ImExRK2CN::mAEx(1,0) = 2./5.;

      // Initialize explicit b factors
      ImExRK2CN::mBEx(1) = 5./6.;
   }

}
}

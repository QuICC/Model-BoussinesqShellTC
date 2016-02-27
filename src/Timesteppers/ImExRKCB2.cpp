/** 
 * @file ImExRKCB2.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 2
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExRKCB2.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Timestep {

   // Scheme requires 3 substeps
   const int ImExRKCB2::STEPS = 3;

   // Scheme order
   const int ImExRKCB2::ORDER = 2;

   // Scheme has embedded lower order scheme?
   const bool ImExRKCB2::HAS_EMBEDDED = true;

   // Name of the scheme
   const std::string ImExRKCB2::NAME = "ImExRKCB2";

   // Set the implicit a factors
   Eigen::Array<MHDFloat,3,3> ImExRKCB2::mAIm = Eigen::Array<MHDFloat,3,3>::Zero();

   // Set the implicit b factors
   Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBIm = Eigen::Array<MHDFloat,3,1>::Zero();

   // Set the explicit a factors
   Eigen::Array<MHDFloat,3,3> ImExRKCB2::mAEx = Eigen::Array<MHDFloat,3,3>::Zero();

   // Set the explicit b factors
   Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBEx = Eigen::Array<MHDFloat,3,1>::Zero();

   // Set the explicit c factors
   Eigen::Array<MHDFloat,3,1> ImExRKCB2::mCEx = Eigen::Array<MHDFloat,3,1>::Zero();

   // Set the implicit embedded b factors
   Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBImErr = Eigen::Array<MHDFloat,3,1>::Zero();

   // Set the explicit embedded b factors
   Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBExErr = Eigen::Array<MHDFloat,3,1>::Zero();

   MHDFloat ImExRKCB2::aIm(const int i, const int j)
   {
      return ImExRKCB2::mAIm(i,j);
   }

   MHDFloat ImExRKCB2::bIm(const int i)
   {
      return ImExRKCB2::mBIm(i);
   }

   MHDFloat ImExRKCB2::aEx(const int i, const int j)
   {
      return ImExRKCB2::mAEx(i,j);
   }

   MHDFloat ImExRKCB2::bEx(const int i)
   {
      return ImExRKCB2::mBEx(i);
   }

   MHDFloat ImExRKCB2::cEx(const int i)
   {
      return ImExRKCB2::mCEx(i);
   }

   MHDFloat ImExRKCB2::bImErr(const int i)
   {
      return ImExRKCB2::mBImErr(i);
   }

   MHDFloat ImExRKCB2::bExErr(const int i)
   {
      return ImExRKCB2::mBExErr(i);
   }

   void ImExRKCB2::init()
   {
      // Initialize implicit a factors
      ImExRKCB2::mAIm(1,1) = 2./5.;
      ImExRKCB2::mAIm(2,1) = 5./6.;
      ImExRKCB2::mAIm(2,2) = 1./6.;

      // Initialize implicit b factors
      ImExRKCB2::mBIm(1) = 5./6.;
      ImExRKCB2::mBIm(2) = 1./6.;

      // Initialize explicit a factors
      ImExRKCB2::mAEx(1,0) = 2./5.;
      ImExRKCB2::mAEx(2,1) = 1.;

      // Initialize explicit b factors
      ImExRKCB2::mBEx(1) = 5./6.;
      ImExRKCB2::mBEx(2) = 1./6.;

      // Initialize explicit c factors
      ImExRKCB2::mCEx(1) = 2./5.;
      ImExRKCB2::mCEx(2) = 1.;
      
      // Initialize implicit embedded b factors
      ImExRKCB2::mBImErr(1) = 4./5.;
      ImExRKCB2::mBImErr(2) = 1./5.;

      // Initialize explicit embedded b factors
      ImExRKCB2::mBExErr(1) = 4./5.;
      ImExRKCB2::mBExErr(2) = 1./5.;
   }

}
}

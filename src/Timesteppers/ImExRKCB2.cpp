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
   const int ImExRKCB2::STEPS = 4;

   // Set the implicit a factors
   const Eigen::Array<MHDFloat,3,3> ImExRKCB2::mAIm = Eigen::Array<MHDFloat,3,3>::Zero();

   // Set the implicit b factors
   const Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBIm = Eigen::Array<MHDFloat,3,1>::Zero();

   // Set the explicit a factors
   const Eigen::Array<MHDFloat,3,3> ImExRKCB2::mAEx = Eigen::Array<MHDFloat,3,3>::Zero();

   // Set the explicit b factors
   const Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBEx = Eigen::Array<MHDFloat,3,1>::Zero();

   // Set the implicit embedded b factors
   const Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBImErr = Eigen::Array<MHDFloat,3,1>::Zero();

   // Set the explicit embedded b factors
   const Eigen::Array<MHDFloat,3,1> ImExRKCB2::mBExErr = Eigen::Array<MHDFloat,3,1>::Zero();

   MHDFloat ImExRKCB2::aIm(const int i, const int j)
   {
      return ImExRKCB2::mAIm(i,j);
   }

   MHDFloat ImExRKCB2::bIm(const int step)
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
      imExRKCB2::aIm(1,1) = 2./5.;
      imExRKCB2::aIm(2,1) = 5./6.;
      imExRKCB2::aIm(2,1) = 1./6.;

      // Initialize implicit b factors
      imExRKCB2::bIm(1) = 5./6.;
      imExRKCB2::bIm(2) = 1./6.;

      // Initialize explicit a factors
      imExRKCB2::aEx(1,0) = 2./5.;
      imExRKCB2::aEx(2,1) = 1.;

      // Initialize explicit b factors
      imExRKCB2::bEx(1) = 5./6.;
      imExRKCB2::bEx(2) = 1./6.;
      
      // Initialize implicit embedded b factors
      imExRKCB2::bImErr(1) = 4./5.;
      imExRKCB2::bImErr(2) = 1./5.;

      // Initialize explicit embedded b factors
      imExRKCB2::bExErr(1) = 4./5.;
      imExRKCB2::bExErr(2) = 1./5.;
   }

}
}

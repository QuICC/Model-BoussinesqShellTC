/** 
 * @file ImExRKCB3b.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3b
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExRKCB3b.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   const int ImExRKCB3b::STEPS = 4;

   // Scheme order
   const int ImExRKCB3b::ORDER = 3;

   // Scheme has embedded lower order scheme?
   const bool ImExRKCB3b::HAS_EMBEDDED = false;

   // Use scheme's embedded lower order scheme?
   bool ImExRKCB3b::USE_EMBEDDED = false;

   // Name of the scheme
   const std::string ImExRKCB3b::NAME = "ImExRKCB3b";

   // Set the implicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3b::mAIm = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the implicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3b::mBIm = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3b::mAEx = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the explicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3b::mBEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit c factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3b::mCEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the implicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3b::mBImErr = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3b::mBExErr = Eigen::Array<MHDFloat,4,1>::Zero();

   MHDFloat ImExRKCB3b::aIm(const int i, const int j)
   {
      return ImExRKCB3b::mAIm(i,j);
   }

   MHDFloat ImExRKCB3b::bIm(const int i)
   {
      return ImExRKCB3b::mBIm(i);
   }

   MHDFloat ImExRKCB3b::aEx(const int i, const int j)
   {
      return ImExRKCB3b::mAEx(i,j);
   }

   MHDFloat ImExRKCB3b::bEx(const int i)
   {
      return ImExRKCB3b::mBEx(i);
   }

   MHDFloat ImExRKCB3b::cEx(const int i)
   {
      return ImExRKCB3b::mCEx(i);
   }

   MHDFloat ImExRKCB3b::bImErr(const int i)
   {
      return ImExRKCB3b::mBImErr(i);
   }

   MHDFloat ImExRKCB3b::bExErr(const int i)
   {
      return ImExRKCB3b::mBExErr(i);
   }

   void ImExRKCB3b::init()
   {
      MHDFloat g = 1./2. + std::sqrt(3.0)/6.0;
      MHDFloat c3 = 1./2. - std::sqrt(3.0)/6.0;
      MHDFloat c4 = 1./2. + std::sqrt(3.0)/6.0;

      // Initialize implicit a factors
      ImExRKCB3b::mAIm(1,1) = g;
      ImExRKCB3b::mAIm(2,1) = -std::sqrt(3.0)/3.0;
      ImExRKCB3b::mAIm(2,2) = g;
      ImExRKCB3b::mAIm(3,3) = g;

      // Initialize implicit b factors
      ImExRKCB3b::mBIm(2) = 1./2.;
      ImExRKCB3b::mBIm(3) = 1./2.;

      // Initialize explicit a factors
      ImExRKCB3b::mAEx(1,0) = g;
      ImExRKCB3b::mAEx(2,1) = c3;
      ImExRKCB3b::mAEx(3,2) = c4;

      // Initialize explicit b factors
      ImExRKCB3b::mBEx(2) = 1./2.;
      ImExRKCB3b::mBEx(3) = 1./2.;

      // Initialize explicit c factors
      ImExRKCB3b::mCEx(1) = g;
      ImExRKCB3b::mCEx(2) = c3;
      ImExRKCB3b::mCEx(3) = c4;
   }

   void ImExRKCB3b::useEmbedded()
   {
      if(ImExRKCB3b::HAS_EMBEDDED)
      {
         ImExRKCB3b::USE_EMBEDDED = true;
      } else
      {
         throw Exception("Tried to activate inexistant embedded scheme!");
      }
   }

}
}

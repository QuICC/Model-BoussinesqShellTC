/** 
 * @file ImExRKCB3e.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3e
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExRKCB3e.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   const int ImExRKCB3e::STEPS = 4;

   // Scheme order
   const int ImExRKCB3e::ORDER = 3;

   // Scheme has embedded lower order scheme?
   const bool ImExRKCB3e::HAS_EMBEDDED = false;

   // Use scheme's embedded lower order scheme?
   bool ImExRKCB3e::USE_EMBEDDED = false;

   // Name of the scheme
   const std::string ImExRKCB3e::NAME = "ImExRKCB3e";

   // Set the implicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3e::mAIm = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the implicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3e::mBIm = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3e::mAEx = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the explicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3e::mBEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit c factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3e::mCEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the implicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3e::mBImErr = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3e::mBExErr = Eigen::Array<MHDFloat,4,1>::Zero();

   MHDFloat ImExRKCB3e::aIm(const int i, const int j)
   {
      return ImExRKCB3e::mAIm(i,j);
   }

   MHDFloat ImExRKCB3e::bIm(const int i)
   {
      return ImExRKCB3e::mBIm(i);
   }

   MHDFloat ImExRKCB3e::aEx(const int i, const int j)
   {
      return ImExRKCB3e::mAEx(i,j);
   }

   MHDFloat ImExRKCB3e::bEx(const int i)
   {
      return ImExRKCB3e::mBEx(i);
   }

   MHDFloat ImExRKCB3e::cEx(const int i)
   {
      return ImExRKCB3e::mCEx(i);
   }

   MHDFloat ImExRKCB3e::bImErr(const int i)
   {
      return ImExRKCB3e::mBImErr(i);
   }

   MHDFloat ImExRKCB3e::bExErr(const int i)
   {
      return ImExRKCB3e::mBExErr(i);
   }

   void ImExRKCB3e::init()
   {
      // Initialize implicit a factors
      ImExRKCB3e::mAIm(1,1) = 1./3.;
      ImExRKCB3e::mAIm(2,1) = 1./2.;
      ImExRKCB3e::mAIm(2,2) = 1./2.;
      ImExRKCB3e::mAIm(3,1) = 3./4.;
      ImExRKCB3e::mAIm(3,2) = -1./4.;
      ImExRKCB3e::mAIm(3,3) = 1./2.;

      // Initialize implicit b factors
      ImExRKCB3e::mBIm(1) = 3./4.;
      ImExRKCB3e::mBIm(2) = -1./4.;
      ImExRKCB3e::mBIm(3) = 1./2.;

      // Initialize explicit a factors
      ImExRKCB3e::mAEx(1,0) = 1./3.;
      ImExRKCB3e::mAEx(2,1) = 1.;
      ImExRKCB3e::mAEx(3,1) = 3./4.;
      ImExRKCB3e::mAEx(3,2) = 1./4.;

      // Initialize explicit b factors
      ImExRKCB3e::mBEx(1) = 3./4.;
      ImExRKCB3e::mBEx(2) = -1./4.;
      ImExRKCB3e::mBEx(3) = 1./2.;

      // Initialize explicit c factors
      ImExRKCB3e::mCEx(1) = 1./3.;
      ImExRKCB3e::mCEx(2) = 1.;
      ImExRKCB3e::mCEx(3) = 1.;
   }

   void ImExRKCB3e::useEmbedded()
   {
      if(ImExRKCB3e::HAS_EMBEDDED)
      {
         ImExRKCB3e::USE_EMBEDDED = true;
      } else
      {
         throw Exception("Tried to activate inexistant embedded scheme!");
      }
   }

}
}

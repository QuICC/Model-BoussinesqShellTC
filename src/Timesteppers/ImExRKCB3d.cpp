/** 
 * @file ImExRKCB3d.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3d
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExRKCB3d.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   const int ImExRKCB3d::STEPS = 4;

   // Scheme order
   const int ImExRKCB3d::ORDER = 3;

   // Scheme has embedded lower order scheme?
   const bool ImExRKCB3d::HAS_EMBEDDED = true;

   // Use scheme's embedded lower order scheme?
   bool ImExRKCB3d::USE_EMBEDDED = false;

   // Name of the scheme
   const std::string ImExRKCB3d::NAME = "ImExRKCB3d";

   // Set the implicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3d::mAIm = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the implicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3d::mBIm = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3d::mAEx = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the explicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3d::mBEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit c factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3d::mCEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the implicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3d::mBImErr = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3d::mBExErr = Eigen::Array<MHDFloat,4,1>::Zero();

   MHDFloat ImExRKCB3d::aIm(const int i, const int j)
   {
      return ImExRKCB3d::mAIm(i,j);
   }

   MHDFloat ImExRKCB3d::bIm(const int i)
   {
      return ImExRKCB3d::mBIm(i);
   }

   MHDFloat ImExRKCB3d::aEx(const int i, const int j)
   {
      return ImExRKCB3d::mAEx(i,j);
   }

   MHDFloat ImExRKCB3d::bEx(const int i)
   {
      return ImExRKCB3d::mBEx(i);
   }

   MHDFloat ImExRKCB3d::cEx(const int i)
   {
      return ImExRKCB3d::mCEx(i);
   }

   MHDFloat ImExRKCB3d::bImErr(const int i)
   {
      return ImExRKCB3d::mBImErr(i);
   }

   MHDFloat ImExRKCB3d::bExErr(const int i)
   {
      return ImExRKCB3d::mBExErr(i);
   }

   void ImExRKCB3d::init()
   {
      // Initialize implicit b factors
      ImExRKCB3d::mBIm(0) = 0.0;
      ImExRKCB3d::mBIm(1) = 355931813527.0/1014712533305.0;
      ImExRKCB3d::mBIm(2) = 709215176366.0/1093407543385.0;
      ImExRKCB3d::mBIm(3) = 755675305.0/1258355728177.0;

      // Initialize implicit a factors
      ImExRKCB3d::mAIm(1,0) = 0.0;
      ImExRKCB3d::mAIm(1,1) = 418884414754.0/469594081263.0;
      ImExRKCB3d::mAIm(2,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAIm(2,1) = -304881946513433262434901.0/718520734375438559540570.0;
      ImExRKCB3d::mAIm(2,2) = 684872032315.0/962089110311.0;
      ImExRKCB3d::mAIm(3,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAIm(3,1) = ImExRKCB3d::mBIm(1);
      ImExRKCB3d::mAIm(3,2) = ImExRKCB3d::mBIm(2);
      ImExRKCB3d::mAIm(3,3) = ImExRKCB3d::mBIm(3);

      // Initialize explicit a factors
      ImExRKCB3d::mAEx(1,0) = 418884414754.0/469594081263.0;
      ImExRKCB3d::mAEx(2,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAEx(2,1) = 214744852859.0/746833870870.0;
      ImExRKCB3d::mAEx(3,0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mAEx(3,1) = ImExRKCB3d::mBIm(1);
      ImExRKCB3d::mAEx(3,2) = 658780719778.0/1014712533305.0;

      // Initialize explicit b factors
      ImExRKCB3d::mBEx(0) = ImExRKCB3d::mBIm(0);
      ImExRKCB3d::mBEx(1) = ImExRKCB3d::mBIm(1);
      ImExRKCB3d::mBEx(2) = ImExRKCB3d::mBIm(2);
      ImExRKCB3d::mBEx(3) = ImExRKCB3d::mBIm(3);

      // Initialize explicit c factors
      ImExRKCB3d::mCEx(1) = ImExRKCB3d::mAEx(1,0);
      ImExRKCB3d::mCEx(2) = ImExRKCB3d::mAEx(2,1);
      ImExRKCB3d::mCEx(3) = 1.0;
      
      // Initialize implicit embedded b factors
      ImExRKCB3d::mBImErr(0) = 0.0;
      ImExRKCB3d::mBImErr(1) = 226763370689.0/646029759300.0;
      ImExRKCB3d::mBImErr(2) = 1496839794860.0/2307829317197.0;
      ImExRKCB3d::mBImErr(3) = 353416193.0/889746336234.0;

      // Initialize explicit embedded b factors
      ImExRKCB3d::mBExErr(0) = 1226988580973.0/2455716303853.0;
      ImExRKCB3d::mBExErr(1) = 0.0;
      ImExRKCB3d::mBExErr(2) = 827818615.0/1665592077861.0;
      ImExRKCB3d::mBExErr(3) = 317137569431.0/634456480332.0;
   }

   void ImExRKCB3d::useEmbedded()
   {
      if(ImExRKCB3d::HAS_EMBEDDED)
      {
         ImExRKCB3d::USE_EMBEDDED = true;
      } else
      {
         throw Exception("Tried to activate inexistant embedded scheme!");
      }
   }

}
}

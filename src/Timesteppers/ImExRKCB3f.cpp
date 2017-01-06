/** 
 * @file ImExRKCB3f.cpp
 * @brief Implementation of an implicit/explicit Runge-Kutta scheme of order 3f
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Timesteppers/ImExRKCB3f.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace QuICC {

namespace Timestep {

   // Scheme requires 3 substeps
   const int ImExRKCB3f::STEPS = 4;

   // Scheme order
   const int ImExRKCB3f::ORDER = 3;

   // Scheme has embedded lower order scheme?
   const bool ImExRKCB3f::HAS_EMBEDDED = true;

   // Use scheme's embedded lower order scheme?
   bool ImExRKCB3f::USE_EMBEDDED = false;

   // Name of the scheme
   const std::string ImExRKCB3f::NAME = "ImExRKCB3f";

   // Set the implicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3f::mAIm = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the implicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3f::mBIm = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit a factors
   Eigen::Array<MHDFloat,4,4> ImExRKCB3f::mAEx = Eigen::Array<MHDFloat,4,4>::Zero();

   // Set the explicit b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3f::mBEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit c factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3f::mCEx = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the implicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3f::mBImErr = Eigen::Array<MHDFloat,4,1>::Zero();

   // Set the explicit embedded b factors
   Eigen::Array<MHDFloat,4,1> ImExRKCB3f::mBExErr = Eigen::Array<MHDFloat,4,1>::Zero();

   MHDFloat ImExRKCB3f::aIm(const int i, const int j)
   {
      return ImExRKCB3f::mAIm(i,j);
   }

   MHDFloat ImExRKCB3f::bIm(const int i)
   {
      return ImExRKCB3f::mBIm(i);
   }

   MHDFloat ImExRKCB3f::aEx(const int i, const int j)
   {
      return ImExRKCB3f::mAEx(i,j);
   }

   MHDFloat ImExRKCB3f::bEx(const int i)
   {
      return ImExRKCB3f::mBEx(i);
   }

   MHDFloat ImExRKCB3f::cEx(const int i)
   {
      return ImExRKCB3f::mCEx(i);
   }

   MHDFloat ImExRKCB3f::bImErr(const int i)
   {
      return ImExRKCB3f::mBImErr(i);
   }

   MHDFloat ImExRKCB3f::bExErr(const int i)
   {
      return ImExRKCB3f::mBExErr(i);
   }

   void ImExRKCB3f::init()
   {
      MHDFloat b1 = -2179897048956./603118880443.;
      MHDFloat b2 = 99189146040./891495457793.;
      MHDFloat b3 = 6064140186914/1415701440113.;
      MHDFloat b4 = 146791865627./668377518349.;

      MHDFloat c2 = 49./50.;
      MHDFloat c3 = 1./25.;

      // Initialize implicit a factors
      ImExRKCB3f::mAIm(1,0) = c2/2.0;
      ImExRKCB3f::mAIm(1,1) = c2/2.0;
      ImExRKCB3f::mAIm(2,0) = -785157464198./1093480182337.;
      ImExRKCB3f::mAIm(2,1) = -30736234873./978681420651.;
      ImExRKCB3f::mAIm(2,2) = 983779726483./1246172347126.;
      ImExRKCB3f::mAIm(3,0) = b1;
      ImExRKCB3f::mAIm(3,1) = b2;
      ImExRKCB3f::mAIm(3,2) = b3;
      ImExRKCB3f::mAIm(3,3) = b4;

      // Initialize implicit b factors
      ImExRKCB3f::mBIm(0) = b1;
      ImExRKCB3f::mBIm(1) = b2;
      ImExRKCB3f::mBIm(2) = b3;
      ImExRKCB3f::mBIm(3) = b4;

      // Initialize explicit a factors
      ImExRKCB3f::mAEx(1,0) = c2;
      ImExRKCB3f::mAEx(2,0) = 13244205847./647648310246.;
      ImExRKCB3f::mAEx(2,1) = 13419997131./686433909488.;
      ImExRKCB3f::mAEx(3,0) = b1;
      ImExRKCB3f::mAEx(3,1) = 231677526244./1085522130027.;
      ImExRKCB3f::mAEx(3,2) = 3007879347537./683461566472.;

      // Initialize explicit b factors
      ImExRKCB3f::mBEx(0) = b1;
      ImExRKCB3f::mBEx(1) = b2;
      ImExRKCB3f::mBEx(2) = b3;
      ImExRKCB3f::mBEx(3) = b4;

      // Initialize explicit c factors
      ImExRKCB3f::mCEx(1) = c2;
      ImExRKCB3f::mCEx(2) = c3;
      ImExRKCB3f::mCEx(3) = 1.;

      // Initialize embedded implicit b factors
      ImExRKCB3f::mBImErr(1) = 337712514207./759004992869.;
      ImExRKCB3f::mBImErr(2) = 311412265155./608745789881.;
      ImExRKCB3f::mBImErr(3) = 52826596233./1214539205236.;

      // Initialize embedded explicit b factors
      ImExRKCB3f::mBExErr(2) = 25./48.;
      ImExRKCB3f::mBExErr(3) = 23./48.;
   }

   void ImExRKCB3f::useEmbedded()
   {
      if(ImExRKCB3f::HAS_EMBEDDED)
      {
         ImExRKCB3f::USE_EMBEDDED = true;
      } else
      {
         throw Exception("Tried to activate inexistant embedded scheme!");
      }
   }

}
}

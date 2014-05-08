/**
 * @file TestTFTDiffusion3D.hpp
 * @brief Implementation of the TFT test equation for 3D diffusion 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTTFTDIFFUSION3D_HPP
#define TESTTFTDIFFUSION3D_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Equations/IScalarEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the TFT test equation for 3D diffusion
    */
   class TestTFTDiffusion3D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pyName     Python script name
          * @param spEqParams  Shared equation parameters
          */
         TestTFTDiffusion3D(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestTFTDiffusion3D();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

      private:
   };

   /// Typedef for a shared TestTFTDiffusion3D
   typedef SharedPtrMacro<TestTFTDiffusion3D> SharedTestTFTDiffusion3D;

}
}

#endif // TESTTFTDIFFUSION3D_HPP

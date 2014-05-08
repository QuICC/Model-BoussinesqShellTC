/**
 * @file TestTFFDiffusion3D.hpp
 * @brief Implementation of the TFF test equation for 3D diffusion 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTTFFDIFFUSION3D_HPP
#define TESTTFFDIFFUSION3D_HPP

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
    * @brief Implementation of the TFF test equation for 3D diffusion
    */
   class TestTFFDiffusion3D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pyName     Python script name
          * @param spEqParams  Shared equation parameters
          */
         TestTFFDiffusion3D(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestTFFDiffusion3D();

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

   /// Typedef for a shared TestTFFDiffusion3D
   typedef SharedPtrMacro<TestTFFDiffusion3D> SharedTestTFFDiffusion3D;

}
}

#endif // TESTTFFDIFFUSION3D_HPP

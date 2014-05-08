/**
 * @file TestTTTDiffusion2D.hpp
 * @brief Implementation of the TTT test equation for 2D diffusion (within 3D model) 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTTTTDIFFUSION2D_HPP
#define TESTTTTDIFFUSION2D_HPP

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
    * @brief Implementation of the TTT test equation for 2D diffusion (within 3D model)
    */
   class TestTTTDiffusion2D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pyName     Python script name
          * @param spEqParams  Shared equation parameters
          */
         TestTTTDiffusion2D(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestTTTDiffusion2D();

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

   /// Typedef for a shared TestTTTDiffusion2D
   typedef SharedPtrMacro<TestTTTDiffusion2D> SharedTestTTTDiffusion2D;

}
}

#endif // TESTTTTDIFFUSION2D_HPP

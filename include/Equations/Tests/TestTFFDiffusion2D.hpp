/**
 * @file TestTFFDiffusion2D.hpp
 * @brief Implementation of the TFF test equation for 2D diffusion (within 3D model) 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTTFFDIFFUSION2D_HPP
#define TESTTFFDIFFUSION2D_HPP

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
    * @brief Implementation of the TFF test equation for 2D diffusion (within 3D model)
    */
   class TestTFFDiffusion2D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         TestTFFDiffusion2D(const std::string& pyNameSharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestTFFDiffusion2D();

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

   /// Typedef for a shared TestTFFDiffusion2D
   typedef SharedPtrMacro<TestTFFDiffusion2D> SharedTestTFFDiffusion2D;

}
}

#endif // TESTTFFDIFFUSION2D_HPP

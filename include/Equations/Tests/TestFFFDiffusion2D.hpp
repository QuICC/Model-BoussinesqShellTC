/**
 * @file TestFFFDiffusion2D.hpp
 * @brief Implementation of the FFF test equation for 2D diffusion (within 3D model) 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TESTFFFDIFFUSION2D_HPP
#define TESTFFFDIFFUSION2D_HPP

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
    * @brief Implementation of the FFF test equation for 2D diffusion (within 3D model)
    */
   class TestFFFDiffusion2D: public IScalarEquation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         TestFFFDiffusion2D(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Destructor
          */
         virtual ~TestFFFDiffusion2D();

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

   /// Typedef for a shared TestFFFDiffusion2D
   typedef SharedPtrMacro<TestFFFDiffusion2D> SharedTestFFFDiffusion2D;

}
}

#endif // TESTFFFDIFFUSION2D_HPP

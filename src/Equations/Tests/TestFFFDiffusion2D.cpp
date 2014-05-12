/** 
 * @file TestFFFDiffusion2D.cpp
 * @brief Source of the implementation of the FFF test equation for 2D diffusion (within 3D model)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tests/TestFFFDiffusion2D.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestFFFDiffusion2D::TestFFFDiffusion2D(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName, spEqParams)
   {
   }

   TestFFFDiffusion2D::~TestFFFDiffusion2D()
   {
   }

   void TestFFFDiffusion2D::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestFFFDiffusion2D::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 0, false, false, false);
   }

   void TestFFFDiffusion2D::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

}
}

/** 
 * @file BoussinesqFPlaneNHBGEStreamfunction.cpp
 * @brief Source of the implementation of the streamfunction equation in the F-plane NHBGE model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlaneNHBGE/Boussinesq/BoussinesqFPlaneNHBGEStreamfunction.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "Enums/FieldIdsTools.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"

namespace QuICC {

namespace Equations {

   BoussinesqFPlaneNHBGEStreamfunction::BoussinesqFPlaneNHBGEStreamfunction(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqFPlaneNHBGEStreamfunction::~BoussinesqFPlaneNHBGEStreamfunction()
   {
   }

   void BoussinesqFPlaneNHBGEStreamfunction::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::PROGNOSTIC, 1, true, false);
   }

   void BoussinesqFPlaneNHBGEStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\nabla^2_{\perp}\psi\f$
      ///
      Physical::StreamAdvection<FieldComponents::Physical::X,FieldComponents::Physical::Y>::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYZ).dom(0).grad(), 1.0);
   }

   void BoussinesqFPlaneNHBGEStreamfunction::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Set solver timing
      this->setSolveTiming(SolveTiming::PROGNOSTIC);

      // Set streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, true));

      // Set non orthogonal vertical vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, true, true));

      //-----------------------------------------------------------------------------------------------------------
      // Restrict components of gradient
      ArrayB   compsVec = ArrayB::Constant(3, true);
      // Don't compute Z component (use fieldIndex get proper index)
      compsVec(fieldIndex(FieldComponents::Physical::Z)) = false;
      std::map<FieldComponents::Spectral::Id,ArrayB>  gradComps;
      gradComps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, compsVec));

      // Update streamfunction gradient requirements
      this->updateFieldRequirements(PhysicalNames::STREAMFUNCTION).updateGradient(gradComps);

      // Update streamfunction gradient requirements
      this->updateFieldRequirements(PhysicalNames::VORTICITYZ).updateGradient(gradComps);

//      //-----------------------------------------------------------------------------------------------------------
//      // Restrict components of 2nd order gradient (only example, actually not used)
//      // Make upper triangular matrix
//      MatrixB   compsTen = MatrixB::Constant(3,3, true);
//      compsTen.triangularView<Eigen::StrictlyLower>().setZero();
//      // Don't compute DxDz derivative (order doesn't matter, but use fieldPairSym)
//      std::pair<int,int> idx = fieldPairSym(FieldComponents::Physical::Z,FieldComponents::Physical::X);
//      compsTen(idx.first, idx.second) = false;
//      // Don't compute DzDy derivative (order doesn't matter, but use fieldPairSym)
//      idx = fieldPairSym(FieldComponents::Physical::Y,FieldComponents::Physical::Z);
//      compsTen(idx.first, idx.second) = false;
//      std::map<FieldComponents::Spectral::Id,MatrixB>  grad2Comps;
//      grad2Comps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, compsTen));
//
//      this->updateFieldRequirements(PhysicalNames::STREAMFUNCTION).updateGradient2(grad2Comps);
   }

}
}

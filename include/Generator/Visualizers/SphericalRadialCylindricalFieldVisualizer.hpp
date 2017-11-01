/**
 * @file SphericalRadialCylindricalFieldVisualizer.hpp
 * @brief Implementation of the spherical vertical component field visualizer 
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

#ifndef SPHERICALRADIALCYLINDRICALFIELDVISUALIZER_HPP
#define SPHERICALRADIALCYLINDRICALFIELDVISUALIZER_HPP

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

namespace QuICC {

namespace Equations {

   /**
    * @brief Implementation of the spherical vertical component field visualizer 
    */
   class SphericalRadialCylindricalFieldVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         SphericalRadialCylindricalFieldVisualizer(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~SphericalRadialCylindricalFieldVisualizer();

         /**
          * @brief Set the vertical field name, the base vector field name and requirements
          */
         void setIdentity(const PhysicalNames::Id vertName, const PhysicalNames::Id fieldName);

         /**
          * @brief Set which component of the ase vector field to use for vertical component (ie. field, gradient, curl)
          */
         void setFieldType(const FieldType::Id type);

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const;

         /**
          * @brief Copy nonlinear calculation back into physical field
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set coupling information
          */
         virtual void setCoupling();

      private:
         /**
          * @brief Storage for output field flag
          */
         FieldType::Id mFieldType;

         /**
          * @brief Storage for output curl flag
          */
         PhysicalNames::Id mFieldName;

         /**
          * @brief Cos(theta)
          */
         Array mCosTheta;

         /**
          * @brief Sin(theta)
          */
         Array mSinTheta;
   };

   /// Typedef for a shared SphericalRadialCylindricalFieldVisualizer
   typedef SharedPtrMacro<SphericalRadialCylindricalFieldVisualizer> SharedSphericalRadialCylindricalFieldVisualizer;

}
}

#endif // SPHERICALRADIALCYLINDRICALFIELDVISUALIZER_HPP

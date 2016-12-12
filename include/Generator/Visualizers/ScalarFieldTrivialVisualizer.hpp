/**
 * @file ScalarFieldTrivialVisualizer.hpp
 * @brief Implementation of the basic scalar field visualizer with trivial solver
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SCALARFIELDTRIVIALVISUALIZER_HPP
#define SCALARFIELDTRIVIALVISUALIZER_HPP

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
    * @brief Implementation of the basic scalar field visualizer with trivial solver
    */
   class ScalarFieldTrivialVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         ScalarFieldTrivialVisualizer(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ScalarFieldTrivialVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient, const bool viewGradient2 = false);

         /**
          * @brief Compute the nonlinear interaction term
          *
          * @param rNLComp Nonlinear term component
          * @param compId  ID of the component (allows for a more general implementation)
          */
         virtual void computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const;

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
         bool mViewField;

         /**
          * @brief Storage for output gradient flag
          */
         bool mViewGradient;

         /**
          * @brief Storage for output 2nd order gradient flag
          */
         bool mViewGradient2;
   };

   /// Typedef for a shared ScalarFieldTrivialVisualizer
   typedef SharedPtrMacro<ScalarFieldTrivialVisualizer> SharedScalarFieldTrivialVisualizer;

}
}

#endif // SCALARFIELDTRIVIALVISUALIZER_HPP

/**
 * @file VectorFieldTrivialVisualizer.hpp
 * @brief Implementation of the basic vector field visualizer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VECTORFIELDTRIVIALVISUALIZER_HPP
#define VECTORFIELDTRIVIALVISUALIZER_HPP

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
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of the basic vector field visualizer
    */
   class VectorFieldTrivialVisualizer: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VectorFieldTrivialVisualizer(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VectorFieldTrivialVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient, const bool viewCurl);

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

         /**
          * @brief Set the nonliner integration components
          */
         virtual void setNLComponents();

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
          * @brief Storage for output curl flag
          */
         bool mViewCurl;
   };

   /// Typedef for a shared VectorFieldTrivialVisualizer
   typedef SharedPtrMacro<VectorFieldTrivialVisualizer> SharedVectorFieldTrivialVisualizer;

}
}

#endif // VECTORFIELDTRIVIALVISUALIZER_HPP

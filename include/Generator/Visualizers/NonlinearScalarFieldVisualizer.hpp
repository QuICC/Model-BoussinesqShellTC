/**
 * @file NonlinearScalarFieldVisualizer.hpp
 * @brief Implementation of a nonlinear scalar field visualizer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef NONLINEARSCALARFIELDVISUALIZER_HPP
#define NONLINEARSCALARFIELDVISUALIZER_HPP

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

   struct NonlinearScalarVisualizerIds
   {
      enum Id {
         CYLINDER_HEAT_ADVECTION = 0,
      };
   };

   /**
    * @brief Implementation of a nonlinear scalar field visualizer
    */
   class NonlinearScalarFieldVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         NonlinearScalarFieldVisualizer(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~NonlinearScalarFieldVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set which fields to output
          */
         void setNonlinearType(const NonlinearScalarVisualizerIds::Id type);

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
          * @brief Type of nonlinear term to compute
          */
         NonlinearScalarVisualizerIds::Id mNonlinearType;
   };

   /// Typedef for a shared NonlinearScalarFieldVisualizer
   typedef SharedPtrMacro<NonlinearScalarFieldVisualizer> SharedNonlinearScalarFieldVisualizer;

}
}

#endif // NONLINEARSCALARFIELDVISUALIZER_HPP

/**
 * @file VectorFieldVisualizer.hpp
 * @brief Implementation of the basic vector field visualizer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VECTORFIELDVISUALIZER_HPP
#define VECTORFIELDVISUALIZER_HPP

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
   class VectorFieldVisualizer: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         VectorFieldVisualizer(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VectorFieldVisualizer();

         /**
          * @brief Set the unknown name and requirements
          */
         void setIdentity(const PhysicalNames::Id name);

         /**
          * @brief Set which fields to output
          */
         void setFields(const bool viewField, const bool viewGradient);

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
   };

   /// Typedef for a shared VectorFieldVisualizer
   typedef SharedPtrMacro<VectorFieldVisualizer> SharedVectorFieldVisualizer;

}
}

#endif // VECTORFIELDVISUALIZER_HPP

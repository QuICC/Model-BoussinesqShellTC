/**
 * @file VelocityStreamVisualizer.hpp
 * @brief Implementation of the streamfunction + axial velocity to velocity field visualizer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VELOCITYSTREAMVISUALIZER_HPP
#define VELOCITYSTREAMVISUALIZER_HPP

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
    * @brief Implementation of the streamfunction + axial velocity to velocity field visualizer
    */
   class VelocityStreamVisualizer: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams Shared equation parameters
          */
         VelocityStreamVisualizer(SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VelocityStreamVisualizer();

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

   /// Typedef for a shared VelocityStreamVisualizer
   typedef SharedPtrMacro<VelocityStreamVisualizer> SharedVelocityStreamVisualizer;

}
}

#endif // VELOCITYSTREAMVISUALIZER_HPP

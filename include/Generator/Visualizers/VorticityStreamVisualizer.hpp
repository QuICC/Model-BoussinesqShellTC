/**
 * @file VorticityStreamVisualizer.hpp
 * @brief Implementation of the streamfunction to vorticity visualizer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VORTICITCYSTREAMVISUALIZER_HPP
#define VORTICITCYSTREAMVISUALIZER_HPP

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
    * @brief Implementation of the streamfunction to vorticity visualizer
    */
   class VorticityStreamVisualizer: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param pyName     Python script name
          * @param spEqParams Shared equation parameters
          */
         VorticityStreamVisualizer(const std::string& pyName, SharedEquationParameters spEqParams);

         /**
          * @brief Simple empty destructor
          */
         virtual ~VorticityStreamVisualizer();

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
          * @brief Set the equation coupling information
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

   /// Typedef for a shared VorticityStreamVisualizer
   typedef SharedPtrMacro<VorticityStreamVisualizer> SharedVorticityStreamVisualizer;

}
}

#endif // VORTICITCYSTREAMVISUALIZER_HPP

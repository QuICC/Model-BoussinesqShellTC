/**
 * @file SphericalTorPolWrapper.hpp
 * @brief Spherical Toroidal/Poloidal decomposition wrapper into velocity field 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERICALTORPOLWRAPPER_HPP
#define SPHERICALTORPOLWRAPPER_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Diagnostics/IVelocityWrapper.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief Spherical Toroidal/Poloidal decomposition wrapper into velocity field 
    */
   class SphericalTorPolWrapper: public IVelocityWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         SphericalTorPolWrapper(const Datatypes::SharedVectorVariableType spTorPol);

         /**
          * @brief Constructor
          */
         ~SphericalTorPolWrapper();

         /**
          * @brief Get first velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& one() const;

         /**
          * @brief Get second velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& two() const;

         /**
          * @brief Get third velocity field component
          */
         virtual const Datatypes::PhysicalScalarType& three() const;

         /**
          * @brief Get Resolution
          */
         virtual const SharedResolution spRes() const;

      protected:

      private:
         /**
          * @brief Shared velocity vector variable
          */
         Datatypes::SharedVectorVariableType mspTorPol;
   };

   /// Typedef for a shared SphericalTorPolWrapper
   typedef SharedPtrMacro<SphericalTorPolWrapper> SharedSphericalTorPolWrapper;
}
}

#endif // SPHERICALTORPOLWRAPPER_HPP

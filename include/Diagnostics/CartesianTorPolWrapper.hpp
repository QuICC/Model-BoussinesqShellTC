/**
 * @file CartesianTorPolWrapper.hpp
 * @brief Cartesian Toroidal/Poloidal decomposition wrapper into velocity field 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIANTORPOLWRAPPER_HPP
#define CARTESIANTORPOLWRAPPER_HPP

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

namespace GeoMHDiSCC {

namespace Diagnostics {

   /**
    * @brief Cartesian Toroidal/Poloidal decomposition wrapper into velocity field 
    */
   class CartesianTorPolWrapper: public IVelocityWrapper
   {
      public:
         /**
          * @brief Constructor
          */
         CartesianTorPolWrapper(const Datatypes::SharedVectorVariableType spTorPol);

         /**
          * @brief Constructor
          */
         ~CartesianTorPolWrapper();

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

   /// Typedef for a shared CartesianTorPolWrapper
   typedef SharedPtrMacro<CartesianTorPolWrapper> SharedCartesianTorPolWrapper;
}
}

#endif // CARTESIANTORPOLWRAPPER_HPP

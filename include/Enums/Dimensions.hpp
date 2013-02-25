/** \file Dimensions.hpp
 *  \brief Definition of some useful enums of the dimensions of the model
 */

#ifndef DIMENSIONS_HPP
#define DIMENSIONS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Dimensions {

   /**
    * @brief Simple struc to label 1D, 2D or 3D simulations
    */
   enum Type {
      /// 1D simulation
      ONED,
      /// 2D simulation
      TWOD,
      /// 3D simulation
      THREED
   };

   /**
    * @brief Simple struct holding the IDs of transform spaces
    */
   struct Transform {
      /**
       * @name Enums for the transform spaces
       */
      enum Id {
         /// First transform space
         TRA1D,
         /// Second transform space
         TRA2D,
         /// Third transform space
         TRA3D
      };
   };

   /**
    * @brief Simple struct holding the IDs of data dimensions
    */
   struct Data {
      /**
       * @name Enums for the data dimensions
       */
      enum Id {
         /// First dimension of data for forward transform
         DATF1D,
         /// First dimension of data for backward transform
         DATB1D,
         /// Second dimension of data
         DAT2D,
         /// Third dimension data
         DAT3D
      };
   };

   /**
    * @brief Simple struct holding IDs for the two spaces
    */
   struct Space {
      /**
       * @name Enums for the dimension spaces
       */
      enum Id {
         /// Spectral space
         SPECTRAL,
         /// Physical space
         PHYSICAL
      };
   };
}
}

#endif // DIMENSIONS_HPP

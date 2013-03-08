/** \file VariableHdf5Tags.hpp
 *  \brief Definitions and names used by the spectral space data files
 */

#ifndef SPECTRALFILEDEFS_HPP
#define SPECTRALFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of HDF5File for a spectral data HDF5 file
    */
   class VariableHdf5Tags
   {
      public:
         /**
          * @brief Truncation group tag name
          */
         static const std::string   TRUNCATION;

         /**
          * @brief Physical truncation group tag name
          */
         static const std::string   TRUNCPHYSICAL;

         /**
          * @brief Spectral truncation group tag name
          */
         static const std::string   TRUNCSPECTRAL;

         /**
          * @brief Truncation dimension tag name
          */
         static const std::string   TRUNCDIM;

         /**
          * @brief Physical parameters part for State file
          */
         static const std::string   PHYSICAL;

      private:
         /**
         * @brief Empty destructor
         */
         VariableHdf5Tags();

         /**
         * @brief Destructor
         */
         ~VariableHdf5Tags();
   };
}
}

#endif // SPECTRALFILEDEFS_HPP

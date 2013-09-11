/** 
 * @file VariableHdf5Tags.hpp
 * @brief Definitions and names used by the spectral space data files
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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

         /**
          * @brief Run parameters part for State file
          */
         static const std::string   RUN;

         /**
          * @brief Time tag for State file
          */
         static const std::string   RUNTIME;

         /**
          * @brief Timestep tag for State file
          */
         static const std::string   RUNSTEP;

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

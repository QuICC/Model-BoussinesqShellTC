/** 
 * @file FScheme.hpp
 * @brief Implementation of the Fourier scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FSCHEME_HPP
#define FSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Splitting.hpp"
#include "Resolutions/Resolution.hpp"
#include "SpatialSchemes/1D/IRegular1DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of Fourier scheme
    */
   class FScheme: public IRegular1DScheme
   {
      public:
         /**
          * @brief Get type string for the scheme
          */
         static std::string type();

         /**
          * @brief Constructor
          *
          * @param dim  Chebyshev truncation
          */
         explicit FScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~FScheme();

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const;

         /**
          * @brief Add the transform setups to resolution
          */
         virtual void addTransformSetups(SharedResolution spRes) const;
         
      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions();

         /**
          * @brief Set transform costs
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setCosts();

         /**
          * @brief Set transform scalings
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setScalings();

         /**
          * @brief Set transform memory footprint
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setMemoryScore();

      private:
         /**
          * @brief Construct setup object for X transform
          */
         Transform::SharedFftSetup  spSetup1D(SharedResolution spRes) const;
   };

}
}

#endif // FSCHEME_HPP

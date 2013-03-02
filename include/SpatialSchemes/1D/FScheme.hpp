/** \file FScheme.hpp
 *  \brief Implementation of the Fourier scheme
 *
 *  \mhdBug Needs test
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
#include "SpatialSchemes/1D/Regular1DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of Fourier scheme
    */
   class FScheme: public Regular1DScheme
   {
      public:
         /**
          * @brief Construct setup object for X transform
          */
         static Transform::SharedFftSetup  spSetup1D(SharedResolution spRes);

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
   };

}

#endif // FSCHEME_HPP

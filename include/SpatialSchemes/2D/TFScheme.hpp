/** \file TFScheme.hpp
 *  \brief Implementation of the Chebyshev(FFT) + Fourier scheme
 *
 *  \mhdBug Needs test
 */

#ifndef TFSCHEME_HPP
#define TFSCHEME_HPP

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
#include "SpatialSchemes/2D/Regular2DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of Chebyshev(FFT) + Fourier scheme
    */
   class TFScheme: public Regular2DScheme
   {
      public:
         /**
          * @brief Construct setup object for first transform
          */
         static Transform::SharedFftSetup  spSetup1D(SharedResolution spRes);

         /**
          * @brief Construct setup object for second transform
          */
         static Transform::SharedFftSetup  spSetup2D(SharedResolution spRes);

         /**
          * @brief Constructor
          *
          * @param dim  Chebyshev truncations
          */
         explicit TFScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~TFScheme();

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
          */
         virtual void setCosts();

         /**
          * @brief Set transform scalings
          */
         virtual void setScalings();

         /**
          * @brief Set transform memory footprint
          */
         virtual void setMemoryScore();

      private:
   };

}

#endif // TFSCHEME_HPP

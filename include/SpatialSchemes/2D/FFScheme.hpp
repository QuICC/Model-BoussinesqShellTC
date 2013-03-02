/** \file FFScheme.hpp
 *  \brief Implementation of the Fourier + Fourier scheme
 *
 *  \mhdBug Needs test
 */

#ifndef FFSCHEME_HPP
#define FFSCHEME_HPP

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
    * @brief Implementation of Fourier + Fourier scheme
    */
   class FFScheme: public Regular2DScheme
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
          * @param dim  Fourier truncations
          */
         explicit FFScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~FFScheme();

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

#endif // FFSCHEME_HPP

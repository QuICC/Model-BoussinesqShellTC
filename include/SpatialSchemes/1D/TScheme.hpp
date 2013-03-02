/** \file TScheme.hpp
 *  \brief Implementation of the Chebyshev(FFT) scheme
 *
 *  \mhdBug Needs test
 */

#ifndef TSCHEME_HPP
#define TSCHEME_HPP

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
    * @brief Implementation of Chebyshev(FFT) scheme
    */
   class TScheme: public Regular1DScheme
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
         explicit TScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~TScheme();

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

#endif // TSCHEME_HPP

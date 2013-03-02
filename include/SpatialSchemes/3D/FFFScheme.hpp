/** \file FFFScheme.hpp
 *  \brief Implementation of the Fourier + Fourier + Fourier scheme
 *
 *  \mhdBug Needs test
 */

#ifndef FFFSCHEME_HPP
#define FFFSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "FastTransforms/FftwTools.hpp"
#include "Enums/Splitting.hpp"
#include "Resolutions/Resolution.hpp"
#include "SpatialSchemes/3D/Regular3DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of Fourier + Fourier + Fourier scheme
    */
   class FFFScheme: public Regular3DScheme
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
          * @brief Construct setup object for third transform
          */
         static Transform::SharedFftSetup  spSetup3D(SharedResolution spRes);

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncations 
          */
         explicit FFFScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~FFFScheme(); 

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

#endif // FFFSCHEME_HPP

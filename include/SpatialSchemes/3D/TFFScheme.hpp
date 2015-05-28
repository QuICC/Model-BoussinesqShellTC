/** 
 * @file TFFScheme.hpp
 * @brief Implementation of the Chebyshev(FFT) + Fourier + Fourier scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TFFSCHEME_HPP
#define TFFSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/FftSelector.hpp"
#include "Enums/Splitting.hpp"
#include "Resolutions/Resolution.hpp"
#include "SpatialSchemes/3D/IRegular3DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of Chebyshev(FFT) + Fourier + Fourier scheme
    */
   class TFFScheme: public IRegular3DScheme
   {
      public:
         /**
          * @brief Interpret the configuration dimensions
          */
         static void interpretConfigDimensions(ArrayI& rDim);

         /**
          * @brief Get type string for the scheme
          */
         static std::string type();

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncations 
          */
         explicit TFFScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~TFFScheme(); 

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const;

         /**
          * @brief Add the transform setups to resolution
          */
         virtual void addTransformSetups(SharedResolution spRes) const;

         /**
          * @brief Add index counter to shared resolution
          */
         virtual void addIndexCounter(SharedResolution spRes);
         
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
         /**
          * @brief Construct setup object for first transform
          */
         Transform::SharedFftSetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::SharedFftSetup  spSetup2D(SharedResolution spRes) const; 

         /**
          * @brief Construct setup object for third transform
          */
         Transform::SharedFftSetup  spSetup3D(SharedResolution spRes) const;

   };

}
}

#endif // TFFSCHEME_HPP

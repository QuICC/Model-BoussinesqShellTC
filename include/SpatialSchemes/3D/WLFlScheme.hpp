/** 
 * @file WLFlScheme.hpp
 * @brief Implementation of the sphere Worland(poly) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral l ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef WLFLSCHEME_HPP
#define WLFLSCHEME_HPP

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
#include "SpatialSchemes/3D/IRegularSHlScheme.hpp"
#include "FastTransforms/FftSetup.hpp"
#include "PolynomialTransforms/PolySetup.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of the sphere Worland(poly) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme with spectral l ordering
    */
   class WLFlScheme: public IRegularSHlScheme
   {
      public:
         /**
          * @brief Get type string for the scheme
          */
         static std::string type();

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncations 
          */
         explicit WLFlScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~WLFlScheme(); 

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
         Transform::SharedPolySetup  spSetup1D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for second transform
          */
         Transform::SharedPolySetup  spSetup2D(SharedResolution spRes) const;

         /**
          * @brief Construct setup object for third transform
          */
         Transform::SharedFftSetup  spSetup3D(SharedResolution spRes) const;
   };

}
}

#endif // WLFLSCHEME_HPP
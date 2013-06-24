/** \file WLFScheme.hpp
 *  \brief Implementation of the spherical Worland + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme
 */

#ifndef WLFSCHEME_HPP
#define WLFSCHEME_HPP

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
#include "SpatialSchemes/3D/IRegularSHScheme.hpp"
#include "FastTransforms/FftSetup.hpp"
#include "PolynomialTransforms/PolySetup.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of the spherical Worland(poly) + Spherical harmonics (Associated Legendre(poly) +  Fourier) scheme
    */
   class WLFScheme: public IRegularSHScheme
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
         explicit WLFScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~WLFScheme(); 

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

#endif // WLFSCHEME_HPP

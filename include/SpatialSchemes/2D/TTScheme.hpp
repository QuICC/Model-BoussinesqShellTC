/** \file TTScheme.hpp
 *  \brief Implementation of the Chebyshev(FFT) + Chebyshev(FFT) scheme
 */

#ifndef TTSCHEME_HPP
#define TTSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/Enums/Splittings.hpp"
#include "Resolutions/Resolution.hpp"
#include "SpatialSchemes/2D/Regular2DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace EPMPhoenix {

   /**
    * @brief Implementation of Chebyshev(FFT) + Chebyshev(FFT) scheme
    */
   class TTScheme: public Regular2DScheme
   {
      public:
         /**
          * @brief Dimensionality of the scheme
          */
         static const int DIMENSIONS;

         /**
          * @brief Construct setup object for first transform
          */
         static SharedFFTSetup  spSetup1D(SharedResolution spRes);

         /**
          * @brief Construct setup object for second transform
          */
         static SharedFFTSetup  spSetup2D(SharedResolution spRes);

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncations
          */
         TTScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~TTScheme() {};

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const;

         /**
          * @brief Get load balancing weights
          */
         virtual Array  loadWeights();

         /**
          * @brief Get memory related score weight
          *
          * @param spRes Resolution information
          */
         virtual double memoryScore(SharedResolution spRes);
         
      protected:
         /**
          * @brief Get size of first truncation
          */
         int nI() const;

         /**
          * @brief Get size of first grid
          */
         int nX() const;

         /**
          * @brief Get size of second truncation
          */
         int nJ() const;

         /**
          * @brief Get size of second grid
          */
         int nY() const;

         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions();

         /**
          * @brief Set transform costs
          */
         virtual void setCosts(const int shift = 0);

         /**
          * @brief Set transform scalings
          */
         virtual void setScalings(const int shift = 0);

         /**
          * @brief Set transform memory footprint
          */
         virtual void setMemory(const int shift = 0);

      private:
   };

   inline int TTScheme::nI() const
   {
      return this->mI + 1;
   }

   inline int TTScheme::nX() const
   {
      return std::ceil((3.0*this->nI())/2.0);
   }

   inline int TTScheme::nJ() const
   {
      return this->mJ + 1;
   }

   inline int TTScheme::nY() const
   {
      return std::ceil((3.0*this->nJ())/2.0);
   }

}

#endif // TTSCHEME_HPP

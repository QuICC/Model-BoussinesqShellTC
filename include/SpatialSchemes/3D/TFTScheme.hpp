/** \file TFTScheme.hpp
 *  \brief Implementation of the Chebyshev(FFT) + Fourier + Chebyshev(FFT) scheme
 */

#ifndef TFTSCHEME_HPP
#define TFTSCHEME_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Transform/FftwTools.hpp"
#include "Base/Enums/Splittings.hpp"
#include "Resolutions/Resolution.hpp"
#include "SpatialSchemes/3D/Regular3DScheme.hpp"
#include "FastTransforms/FftSetup.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of Chebyshev(FFT) + Fourier + Chebyshev(FFT) scheme
    */
   class TFTScheme: public Regular3DScheme
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
          * @brief Construct setup object for third transform
          */
         static SharedFFTSetup  spSetup3D(SharedResolution spRes);

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncations 
          */
         TFTScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~TFTScheme(); 

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
          * @brief Get size of third truncation
          */
         int nK() const;

         /**
          * @brief Get size of third grid
          */
         int nZ() const;

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

   inline int TFTScheme::nI() const
   {
      return this->mI + 1;
   }

   inline int TFTScheme::nX() const
   {
      // Get standard dealiased FFT size
      int nFFT = FftwTools::dealias(this->nI());

      // Check for optimised FFT sizes
      nFFT = FftwTools::optimiseFft(nFFT);

      return nFFT;
   }

   inline int TFTScheme::nJ() const
   {
      return this->mJ + 1;
   }

   inline int TFTScheme::nY() const
   {
      // Get standard dealiased size for mixed FFT
      int nFFT = FftwTools::dealiasMixed(this->nJ());

      // Check for optimised FFT sizes
      nFFT = FftwTools::optimiseFft(nFFT);

      return nFFT;
   }

   inline int TFTScheme::nK() const
   {
      return this->mK + 1;
   }

   inline int TFTScheme::nZ() const
   {
      // Get standard dealiased FFT size
      int nFFT = FftwTools::dealias(this->nK());

      // Check for optimised FFT sizes
      nFFT = FftwTools::optimiseFft(nFFT);

      return nFFT;
   }

}

#endif // TFTSCHEME_HPP

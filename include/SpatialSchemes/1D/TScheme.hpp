/** \file TScheme.hpp
 *  \brief Implementation of the Chebyshev scheme
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
#include "Base/General/Typedefs.hpp"
#include "Base/Enums/Splittings.hpp"
#include "Base/Resolutions/Resolution.hpp"
#include "Base/SpatialSchemes/SpatialScheme.hpp"
#include "Base/Transforms/FFT/FFTSetup.hpp"

namespace EPMPhoenix {

   /**
    * @brief Implementation of Chebyshev scheme
    */
   class TScheme: public SpatialScheme
   {
      public:
         /**
          * @brief Dimensionality of the scheme
          */
         static const int DIMENSIONS;

         /**
          * @brief Construct setup object for X transform
          */
         static SharedFFTSetup  spSetup1D(SharedResolution spRes);

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncation
          * @param shift   Shift of the dimensions
          */
         TScheme(const ArrayI& dim, const int shift = 0);

         /**
          * @brief Destructor
          */
         virtual ~TScheme() {};

         /**
          * @brief Create indexes for a possibly restricted set
          *
          * @param dim     Dimension for which to compute indexes
          * @param fwd1D   Storage for forward indexes of first dimension
          * @param bwd1D   Storage for backward indexes of first dimension
          * @param idx2D   Storage for the indexes of second dimension
          * @param idx3D   Storage for forward indexes of third dimension
          * @param id      ID of the bin
          * @param bins    Total number of bins (useful to build efficient pairs)
          * @param n0      Starting index of restricted set
          * @param nN      Length of restricted set
          * @param flag    Flag to specify location of splitting
          */
         virtual void fillIndexes(const int dim, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id = ArrayI(), const ArrayI& bins = ArrayI(), const ArrayI& n0 = ArrayI(), const ArrayI& nN = ArrayI(), Splittings::Locations::Id flag = Splittings::Locations::NONE);

         /**
          * @brief Get total of splittable indexes 
          *
          * @param dim     Dimension for which to compute indexes
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const int dim, Splittings::Locations::Id flag);

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const;

         /**
          * @brief Get load balancing weights
          */
         virtual Array loadWeights();

         /**
          * @brief Get memory related score weight
          *
          * @param spRes Resolution information
          */
         virtual double memoryScore(SharedResolution spRes);
         
      protected:
         /**
          * @brief Get size of X truncation
          */
         int nI() const;

         /**
          * @brief Get size of X grid
          */
         int nX() const;

         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions();

         /**
          * @brief Set transform costs
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setCosts(const int shift = 0);

         /**
          * @brief Set transform scalings
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setScalings(const int shift = 0);

         /**
          * @brief Set transform memory footprint
          *
          * @param shift   Shift of the dimensions
          */
         virtual void setMemory(const int shift = 0);

      private:
         /**
          * @brief Radial truncation
          */
         int   mI;
   };

   inline int TScheme::nI() const
   {
      return this->mI + 1;
   }

   inline int TScheme::nX() const
   {
      return std::ceil((3.0*this->nI())/2.0);
   }

}

#endif // TSCHEME_HPP

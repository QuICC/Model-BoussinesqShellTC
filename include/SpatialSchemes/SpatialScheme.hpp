/** \file SpatialScheme.hpp
 *  \brief Implementation of the basic components of the spatial scheme
 */

#ifndef SPATIALSCHEME_HPP
#define SPATIALSCHEME_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <assert.h>
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Splitting.hpp"
#include "SpatialSchemes/SchemeBase.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the basic components of the spatial scheme
    */
   class SpatialScheme : public SchemeBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         explicit SpatialScheme(const int dims);

         /**
          * @brief Destructor
          */
         virtual ~SpatialScheme();

         /**
          * @brief Initialise the scheme
          */
         void init();

         /**
          * @brief Create indexes for a possibly restricted set (to simplify implementation is is defined for 3D cases)
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
         virtual void fillIndexes(const int dim, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id = ArrayI(), const ArrayI& bins = ArrayI(), const ArrayI& n0 = ArrayI(), const ArrayI& nN = ArrayI(), Splitting::Locations::Id flag = Splitting::Locations::NONE) = 0;

         /**
          * @brief Get total of splittable indexes 
          *
          * @param dim     Dimension for which to compute indexes
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const int dim, Splitting::Locations::Id flag) = 0;

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const = 0;
         
      protected:
         /**
          * @brief Initialise the domain dimensions
          */
         virtual void setDimensions() = 0;

         /**
          * @brief Get forward dimension of given transform
          *
          * @param i Transform index
          */
         int dimFwd(const int i) const;

         /**
          * @brief Get backward dimension of given transform
          *
          * @param i Transform index
          */
         int dimBwd(const int i) const;

         /**
          * @brief Get second dimension of given transform
          *
          * @param i Transform index
          */
         int dim2D(const int i) const;

         /**
          * @brief Get third dimension of given transform
          *
          * @param i Transform index
          */
         int dim3D(const int i) const;

         /**
          * @brief Get dimension of the domain
          */
         int dims() const;

         /**
          * @brief Full dimensions of the domain
          */
         std::vector<ArrayI>   mDimensions;

      private:
         /**
          * @brief Dimension of the domain
          */
         int   mDims;
   };

   inline int SpatialScheme::dimFwd(const int i) const
   {
      // Assert for domain dimensions
      assert(i < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(i).size() > 0);

      return this->mDimensions.at(i)(0);
   }

   inline int SpatialScheme::dimBwd(const int i) const
   {
      // Assert for domain dimensions
      assert(i < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(i).size() > 1);

      return this->mDimensions.at(i)(1);
   }

   inline int SpatialScheme::dim2D(const int i) const
   {
      // Assert for domain dimensions
      assert(i < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(i).size() > 2);

      return this->mDimensions.at(i)(2);
   }

   inline int SpatialScheme::dim3D(const int i) const
   {
      // Assert for domain dimensions
      assert(i < this->mDims);

      // Assert for dimension size
      assert(this->mDimensions.at(i).size() > 3);

      return this->mDimensions.at(i)(3);
   }

   inline int SpatialScheme::dims() const
   {
      return this->mDims;
   }

   /// Typedef for a shared pointer to a SpatialScheme object
   typedef SharedPtrMacro<SpatialScheme>   SharedSpatialScheme;
}


#endif // SPATIALSCHEME_HPP

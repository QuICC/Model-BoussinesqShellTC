/** 
 * @file ISpatialScheme.hpp
 * @brief Implementation of the basic components of the spatial scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ISPATIALSCHEME_HPP
#define ISPATIALSCHEME_HPP

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
#include "SpatialSchemes/ISchemeCosts.hpp"
#include "LoadSplitter/Algorithms/SplittingDescription.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of the basic components of the spatial scheme
    */
   class ISpatialScheme : public ISchemeCosts
   {
      public:
         /**
          * @brief Interpret the configuration dimensions
          */
         static void interpretConfigDimensions(ArrayI& dim);

         /**
          * @brief Tune the shared resolution used by simulation
          */
         static void tuneResolution(SharedResolution spRes, const Parallel::SplittingDescription& descr);

         /**
          * @brief Constructor
          *
          * @param dims Dimension of the domain
          */
         explicit ISpatialScheme(const int dims);

         /**
          * @brief Destructor
          */
         virtual ~ISpatialScheme();

         /**
          * @brief Initialise the scheme
          */
         void init();

         /**
          * @brief Create indexes for a possibly restricted set (to simplify implementation is is defined for 3D cases)
          *
          * @param transId Transform ID
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
         virtual int fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id = ArrayI(), const ArrayI& bins = ArrayI(), const ArrayI& n0 = ArrayI(), const ArrayI& nN = ArrayI(), const Splitting::Locations::Id flag = Splitting::Locations::NONE) = 0;

         /**
          * @brief Get total of splittable indexes 
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag) = 0;

         /**
          * @brief Scheme specific splitting restrictions
          */
         virtual bool applicable() const = 0;

         /**
          * @brief Add scheme specific transform setups to resolution
          */
         virtual void addTransformSetups(SharedResolution spRes) const = 0;

         /**
          * @brief Get the simulation wide spectral array dimensions (can be different from spectral resolution)
          */
         const ArrayI& getTransformSpace() const;

         /**
          * @brief Add index counter to shared resolution
          */
         virtual void addIndexCounter(SharedResolution spRes);
         
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
         int dim(const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId) const;

         /**
          * @brief Get forward dimension of given transform
          *
          * @param i Transform index
          */
         void setDimension(int d, const Dimensions::Transform::Id transId, const Dimensions::Data::Id dataId);

         /**
          * @brief Set transform dimension
          */
         void setTransformSpace(const ArrayI& dim);

         /**
          * @brief Get dimension of the domain
          */
         int dims() const;

         /**
          * @brief Tune resolution with MPI related conditions
          */
         static void tuneMpiResolution(const Parallel::SplittingDescription& descr);

      private:
         /**
          * @brief Dimension of the domain
          */
         int   mDims;

         /**
          * @brief Array size of the spectral arrays (can be different from spectral resolution)
          */
         ArrayI mTransformSpace;

         /**
          * @brief Full dimensions of the domain
          */
         std::vector<ArrayI>   mDimensions;
   };

   /// Typedef for a shared pointer to a ISpatialScheme object
   typedef SharedPtrMacro<ISpatialScheme>   SharedISpatialScheme;
}
}

#endif // ISPATIALSCHEME_HPP

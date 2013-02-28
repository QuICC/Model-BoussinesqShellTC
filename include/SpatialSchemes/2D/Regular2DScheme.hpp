/** \file Regular2DScheme.hpp
 *  \brief Implementation of a generic regular 2D scheme
 *
 *  \mhdBug Needs test
 */

#ifndef REGULAR2DSCHEME_HPP
#define REGULAR2DSCHEME_HPP

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
#include "SpatialSchemes/SpatialScheme.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of a generic regular 2D scheme
    */
   class Regular2DScheme: public SpatialScheme
   {
      public:
         /**
          * @brief Dimensionality of the scheme
          */
         static const int DIMENSIONS;

         /**
          * @brief Constructor
          *
          * @param dim     Chebyshev truncations
          */
         explicit Regular2DScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~Regular2DScheme();

         /**
          * @brief Create indexes for a possibly restricted set
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
         virtual void fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id = ArrayI(), const ArrayI& bins = ArrayI(), const ArrayI& n0 = ArrayI(), const ArrayI& nN = ArrayI(), Splitting::Locations::Id flag = Splitting::Locations::NONE);

         /**
          * @brief Get total of splittable indexes 
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag);
         
      protected:
         /**
          * @brief First truncation
          */
         int   mI;

         /**
          * @brief Second truncation
          */
         int   mJ;

      private:
   };

}

#endif // REGULAR2DSCHEME_HPP

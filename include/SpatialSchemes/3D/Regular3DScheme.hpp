/** \file Regular3DScheme.hpp
 *  \brief Implementation of a generic regular 3D scheme
 */

#ifndef REGULAR3DSCHEME_HPP
#define REGULAR3DSCHEME_HPP

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
#include "SpatialSchemes/SpatialScheme.hpp"

namespace EPMPhoenix {

   /**
    * @brief Implementation of a generic regular 3D scheme
    */
   class Regular3DScheme: public SpatialScheme
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
         Regular3DScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~Regular3DScheme() {};

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
         
      protected:
         /**
          * @brief First regular truncation
          */
         int   mI;

         /**
          * @brief Second regular truncation
          */
         int   mJ;

         /**
          * @brief third regular truncation
          */
         int   mK;

      private:
   };

}

#endif // REGULAR3DSCHEME_HPP

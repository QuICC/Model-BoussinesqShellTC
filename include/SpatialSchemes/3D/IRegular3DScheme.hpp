/** 
 * @file IRegular3DScheme.hpp
 * @brief Implementation of a generic regular 3D scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef IREGULAR3DSCHEME_HPP
#define IREGULAR3DSCHEME_HPP

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
#include "SpatialSchemes/ISpatialScheme.hpp"

namespace GeoMHDiSCC {

namespace Schemes {

   /**
    * @brief Implementation of a generic regular 3D scheme
    */
   class IRegular3DScheme: public ISpatialScheme
   {
      public:
         /**
          * @brief Dimensionality of the scheme
          */
         static const int DIMENSIONS;

         /**
          * @brief Data regularity flag
          */
         static bool isRegular();

         /**
          * @brief Constructor
          *
          * @param dim Dimension truncations 
          */
         explicit IRegular3DScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~IRegular3DScheme();

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
}

#endif // IREGULAR3DSCHEME_HPP

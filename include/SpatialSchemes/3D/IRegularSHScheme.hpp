/** 
 * @file IRegularSHScheme.hpp
 * @brief Implementation of the Regular basis + Spherical Harmonics scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IREGULARSHSCHEME_HPP
#define IREGULARSHSCHEME_HPP

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
    * @brief Implementation of Regular basis + Spherical Harmonics scheme
    */
   class IRegularSHScheme: public ISpatialScheme
   {
      public:
         /**
          * @brief Tune the resolution
          */
         static void tuneResolution(SharedResolution spRes);

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
          * @param dim 0: radial, 1: latitudinal, 2: longitudinal
          */
         explicit IRegularSHScheme(const ArrayI& dim);

         /**
          * @brief Destructor
          */
         virtual ~IRegularSHScheme();

         /**
          * @brief Create indexes for a possibly restricted set
          *
          * \mhdTodo Partial splitting needs to be optimized
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
          * @brief Regular truncation
          */
         int   mI;

         /**
          * @brief Spherical harmonic degree
          */
         int   mL;

         /**
          * @brief Spherical harmonic order
          */
         int   mM;

      private:
         /**
          * @brief Build the ML map of the spherical harmonics
          */
         void buildMLMap(std::multimap<int,int>& harmonics, const int id, const int bins);
   };
}
}

#endif // IREGULARSHSCHEME_HPP

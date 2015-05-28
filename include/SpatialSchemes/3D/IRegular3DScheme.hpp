/** 
 * @file IRegular3DScheme.hpp
 * @brief Implementation of a generic regular 3D scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
         virtual int fillIndexes(const Dimensions::Transform::Id transId, std::vector<ArrayI>& fwd1D, std::vector<ArrayI>& bwd1D, std::vector<ArrayI>& idx2D, ArrayI& idx3D, const ArrayI& id = ArrayI(), const ArrayI& bins = ArrayI(), const ArrayI& n0 = ArrayI(), const ArrayI& nN = ArrayI(), const Splitting::Locations::Id flag = Splitting::Locations::NONE);

         /**
          * @brief Get total of splittable indexes 
          *
          * @param transId Transform ID
          * @param flag    Flag to specify location of splitting
          */
         virtual int splittableTotal(const Dimensions::Transform::Id transId, Splitting::Locations::Id flag);
         
      protected:
         /**
          * @brief Compute mode distribution for serial algorithm
          *
          * @param modes   Map of all modes
          */
         void splitSerial(std::multimap<int,int>& modes, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for single splitting on 1D algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitSingle1D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for single splitting on 2D algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitSingle2D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for coupled 2D algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitCoupled2D(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

         /**
          * @brief Compute mode distribution for tubular 2D decomposition algorithm
          *
          * @param modes   Map of all modes
          * @param n0      Starting mode
          * @param nN      Number of modes
          */
         void splitTubular(std::multimap<int,int>& modes, const ArrayI& n0, const ArrayI& nN, const Dimensions::Transform::Id transId);

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

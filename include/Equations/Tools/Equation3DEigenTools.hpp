/**
 * @file Equation3DEigenTools.hpp
 * @brief Implementation of some tools for schemes with three eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATION3DEIGENTOOLS_HPP
#define EQUATION3DEIGENTOOLS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Equations/IEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of some tools for schemes with three eigen direction
    */
   class Equation3DEigenTools
   {
      public:
         /**
          * @brief Number of eigen dimensions
          */
         static const int EIGEN_DIMS = 3;

         /**
          * @brief Create setup for matrices with minial coupling
          *
          * @param spRes   Shared resolution
          * @param nMat    Number of matrices
          * @param blocks  Size of the blocks
          * @param cols    Number of right-hand sides
          */
         static void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols);

         /**
          * @brief General implementation of linear row
          */
         template <typename TEquation> static DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of time row
          */
         template <typename TEquation> static DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of boundary row
          */
         template <typename TEquation> static DecoupledZSparse boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of the boundary block
          */
         static void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D);
   };

   template <typename TEquation> DecoupledZSparse Equation3DEigenTools::linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      throw Exception("Not yet implemented!");
   }

   template <typename TEquation> DecoupledZSparse Equation3DEigenTools::timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      throw Exception("Not yet implemented!");
   }

   template <typename TEquation> DecoupledZSparse Equation3DEigenTools::boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      throw Exception("Not yet implemented!");
   }

}
}

#endif // EQUATION3DEIGENTOOLS_HPP

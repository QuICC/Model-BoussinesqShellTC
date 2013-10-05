/**
 * @file EquationEigen3DTools.hpp
 * @brief Implementation of some tools for schemes with three eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGEN3DTOOLS_HPP
#define EQUATIONEIGEN3DTOOLS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include<vector>

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/BoundaryMethodSelector.hpp"
#include "Equations/IEquation.hpp"
#include "Equations/Tools/EquationEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Tools for equations with three eigen directions
 */
namespace Eigen3D {

   /// Flag to specify index independent boundary conditions
   const Boundary::BCIndex INDEPENDENT;

   /**
    * @brief Set eigen values
    */
   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx);

   /**
    * @brief Create setup for matrices with minial coupling
    *
    * @param spRes   Shared resolution
    * @param nMat    Number of matrices
    * @param blocks  Size of the blocks
    * @param cols    Number of right-hand sides
    */
   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols);

   /**
    * @brief General implementation of linear row
    */
   template <typename TEquation> DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief General implementation of time row
    */
   template <typename TEquation> DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief General implementation of boundary row
    */
   template <typename TEquation> void boundaryRow(TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   template <typename TEquation> void storeBoundaryCondition(TEquation& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& coeffs, const std::vector<Boundary::BCIndex>& bcIdx);
//
//   /**
//    * @brief General implementation of the boundary block
//    */
//   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D);


//
// Implementation follows
//

   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx)
   {
      throw Exception("Not yet implemented!");
      std::vector<MHDFloat> eigs;

      // Fill eigs somehow

      return eigs;
   }

   template <typename TEquation> DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      return EigenTools::makeLinearRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      return EigenTools::makeTimeRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> void boundaryRow(TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      EigenTools::makeBoundaryRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> void storeBoundaryCondition(IEquation& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& coeffs, const std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(coeffs.size() == bcIdx.size());
      assert(coeffs.size() == 0);
   }

}
}
}

#endif // EQUATIONEIGEN3DTOOLS_HPP

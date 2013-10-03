/**
 * @file Equation2DEigenTools.hpp
 * @brief Implementation of some tools for schemes with two eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATION2DEIGENTOOLS_HPP
#define EQUATION2DEIGENTOOLS_HPP

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
#include "Equations/Tools/EquationEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of some tools for schemes with two eigen direction
    */
   class Equation2DEigenTools
   {
      public:
         /**
          * @brief Set eigen values
          */
         template <typename TEquation> static std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx);

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

   template <typename TEquation> std::vector<MHDFloat> Equation2DEigenTools::getEigs(const TEquation& eq, const int matIdx)
   {
      std::vector<MHDFloat> eigs;

      // Get mode indexes
      ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // k2D_
      eigs.push_back(eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0)));
      // k3D_
      eigs_push_back(eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1)));

      return eigs;
   }

   template <typename TEquation> DecoupledZSparse Equation2DEigenTools::linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = Equation2DEigenTools::getEigs(eq, matIdx);

      return EquationEigenTools::makeLinearRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> DecoupledZSparse Equation2DEigenTools::timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = Equation2DEigenTools::getEigs(eq, matIdx);

      return EquationEigenTools::makeTimeRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> DecoupledZSparse Equation2DEigenTools::boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = Equation2DEigenTools::getEigs(eq, matIdx);

      return EquationEigenTools::makeBoundaryRow(eq, compId, matIdx, eigs);
   }

}
}

#endif // EQUATION2DEIGENTOOLS_HPP

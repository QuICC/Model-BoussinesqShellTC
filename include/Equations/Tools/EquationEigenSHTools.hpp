/**
 * @file EquationEigenSHTools.hpp
 * @brief Implementation of some tools for schemes with spherical harmonics expansions 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGENSHTOOLS_HPP
#define EQUATIONEIGENSHTOOLS_HPP

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
#include "Base/KroneckerTypedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Equations/IEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Tools for equations with spherical harmonic expansions
 */
namespace EigenSH {

   /**
    * @brief Set eigen values
    */
   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx);

   /**
    * @brief Compute the number of matrix operators for field coupling
    *
    * @param spRes   Shared resolution
    */
   int fieldCouplingNMat(const SharedResolution spRes);

//
// Implementation follows
//

   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx)
   {
      std::vector<MHDFloat> eigs;

      eigs.push_back(static_cast<MHDFloat>(eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx)));

      return eigs;
   }

}
}
}

#endif // EQUATIONEIGENSHTOOLS_HPP

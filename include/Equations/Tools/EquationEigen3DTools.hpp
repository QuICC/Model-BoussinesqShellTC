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

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "Equations/IEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Tools for equations with three eigen directions
 */
namespace Eigen3D {

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
      throw Exception("Not yet implemented!");
      std::vector<MHDFloat> eigs;

      // Fill eigs somehow

      return eigs;
   }

}
}
}

#endif // EQUATIONEIGEN3DTOOLS_HPP

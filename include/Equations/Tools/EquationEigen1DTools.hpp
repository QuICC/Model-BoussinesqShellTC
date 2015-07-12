/**
 * @file EquationEigen1DTools.hpp
 * @brief Implementation of some tools for schemes with a single eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGEN1DTOOLS_HPP
#define EQUATIONEIGEN1DTOOLS_HPP

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

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Tools for equations with a single eigen direction
 */
namespace Eigen1D {

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

   /**
    * @brief Interpret Tau resolution provided by python code
    */
   void interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution spRes);

   /**
    * @brief Interpret Galerkin resolution provided by python code
    */
   void interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution spRes);

   /**
    * @brief Interpret number of RHS provided by python code
    */
   void interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution spRes);

   /**
    * @brief Interpret system size provided by python code
    */
   void interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution spRes);

//
// Implementation follows
//

   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx)
   {
      std::vector<MHDFloat> eigs;

      // Get wave number rescale to box size
      eigs.push_back(eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DATND>(matIdx)));

      return eigs;
   }

}
}
}

#endif // EQUATIONEIGEN1DTOOLS_HPP

/** 
 * @file EquationEigenSHmTools.cpp
 * @brief Source of the tools for schemes with spherical harmonic expansions with m spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>

// Class include
//
#include "Equations/Tools/EquationEigenSHmTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace EigenSHm {

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

   void interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution spRes)
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #ifdef GEOMHDISCC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rTauNs(m) = tauSize*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rTauNs(m) = tauSize*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //GEOMHDISCC_MPISPSOLVE

      }
   }

   void interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution spRes)
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #ifdef GEOMHDISCC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rGalerkinNs(m) = galerkinSize*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rGalerkinNs(m) = galerkinSize*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //GEOMHDISCC_MPISPSOLVE
      }
   }

   void interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution spRes)
   {
      rRhsCols.setConstant(rhsSize);
   }

   void interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution spRes)
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #ifdef GEOMHDISCC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rSystemNs(m) = systemSize*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rSystemNs(m) = systemSize*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //GEOMHDISCC_MPISPSOLVE
      }
   }

}
}
}

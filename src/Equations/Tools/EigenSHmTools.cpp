/** 
 * @file EigenSHmTools.cpp
 * @brief Source of the tools for schemes with spherical harmonic expansions with m spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Tools/EigenSHmTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   EigenSHmTools::EigenSHmTools()
   {
   }

   EigenSHmTools::~EigenSHmTools()
   {
   }

   std::vector<MHDFloat> EigenSHmTools::identifyEigs(const SharedResolution& spRes, const int matIdx) const
   {
      std::vector<MHDFloat> eigs;

      eigs.push_back(static_cast<MHDFloat>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx)));

      return eigs;
   }

   int EigenSHmTools::computeNMat(const SharedResolution& spRes) const
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

   void EigenSHmTools::interpretTauN(ArrayI& rTauNs, const int tauSize, const SharedResolution& spRes) const
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rTauNs(m) = tauSize*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rTauNs(m) = tauSize*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //GEOMHDISCC_MPISPSOLVE

      }
   }

   void EigenSHmTools::interpretGalerkinN(ArrayI& rGalerkinNs, const int galerkinSize, const SharedResolution& spRes) const
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rGalerkinNs(m) = galerkinSize*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rGalerkinNs(m) = galerkinSize*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //GEOMHDISCC_MPISPSOLVE
      }
   }

   void EigenSHmTools::interpretRhsN(ArrayI& rRhsCols, const int rhsSize, const SharedResolution& spRes) const
   {
      rRhsCols.setConstant(rhsSize);
   }

   void EigenSHmTools::interpretSystemN(ArrayI& rSystemNs, const int systemSize, const SharedResolution& spRes) const
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined GEOMHDISCC_MPI && defined GEOMHDISCC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rSystemNs(m) = systemSize*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rSystemNs(m) = systemSize*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //GEOMHDISCC_MPISPSOLVE
      }
   }

}
}

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

namespace QuICC {

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

      eigs.push_back(static_cast<MHDFloat>(spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx)));

      return eigs;
   }

   int EigenSHmTools::computeNMat(const SharedResolution& spRes) const
   {
      return spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
   }

   void EigenSHmTools::interpretTauN(ArrayI& rTauNs, const SharedResolution& spRes) const
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rTauNs(m) = rTauNs(m)*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rTauNs(m) = rTauNs(m)*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //QUICC_MPISPSOLVE

      }
   }

   void EigenSHmTools::interpretGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rGalerkinNs(m) = rGalerkinNs(m)*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rGalerkinNs(m) = rGalerkinNs(m)*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //QUICC_MPISPSOLVE
      }
   }

   void EigenSHmTools::interpretRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const
   {
      // Python setup is sufficient
   }

   void EigenSHmTools::interpretSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const
   {
      for(int m = 0; m < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); m++)
      {
         #if defined QUICC_MPI && defined QUICC_MPISPSOLVE
            int m_ = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m);
            rSystemNs(m) = rSystemNs(m)*(spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL)-m_);
         #else
            rSystemNs(m) = rSystemNs(m)*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(m);
         #endif //QUICC_MPISPSOLVE
      }
   }

}
}

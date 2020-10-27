/**
 * @file SphereNusseltWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>
#include <stdexcept>

// External includes
//

// Class include
//
#include "IoVariable/SphereNusseltWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "Base/MathConstants.hpp"
#include "IoVariable/NusseltTags.hpp"
#include "PolynomialTransforms/WorlandPolynomial.hpp"

namespace QuICC {

namespace IoVariable {

   SphereNusseltWriter::SphereNusseltWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + NusseltTags::BASENAME, NusseltTags::EXTENSION, prefix + NusseltTags::HEADER, type, NusseltTags::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mHasMOrdering(false), mNusselt(std::numeric_limits<MHDFloat>::quiet_NaN()), mTb(std::numeric_limits<MHDFloat>::quiet_NaN()), mOrigin(0,0)
   {
   }

   SphereNusseltWriter::~SphereNusseltWriter()
   {
   }

   void SphereNusseltWriter::init()
   {
      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         this->mHasMOrdering = true;
      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         this->mHasMOrdering = false;
      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      int m0, l0;
      if(this->mHasMOrdering)
      {
         m0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         l0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0);
      } else
      {
         l0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         m0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0);
      }

      // Background state
      this->mTb = 0.5;

      // Look for l = 0, m = 0 mode
      if(m0 == 0 && l0 == 0)
      {
         internal::Array grid = internal::Array::Zero(1);
         int nN = this->res().sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL);
         Matrix poly(grid.size(), nN);
         internal::Matrix  ipoly(grid.size(), nN);
         Polynomial::WorlandPolynomial::Wnl(poly, ipoly, 0, grid);
         this->mOrigin = poly.transpose();
         this->mOrigin /= std::sqrt(4.0*Math::PI);
      } else
      {
         this->mOrigin.resize(0,0);
      }

      IVariableAsciiWriter::init();
   }

   void SphereNusseltWriter::write()
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      if(this->mOrigin.size() > 0)
      {
         this->mNusselt = this->mTb/(this->mTb + (this->mOrigin.transpose()*sRange.first->second->dom(0).total().profile(0,0).real())(0,0));
      } else
      {
         this->mNusselt = 0.0;
      }

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mNusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mNusselt << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mNusselt))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Sphere Nusselt is NaN!");
      }
   }

}
}

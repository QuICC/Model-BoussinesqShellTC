/** \file CylindricalChebyshevFftwTransform.cpp
 *  \brief Source of the implementation of the Chebyshev FFTW transform for a cylindrical radius
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "FastTransforms/CylindricalChebyshevFftwTransform.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/MathConstants.hpp"
#include "FastTransforms/FftwLibrary.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   Array CylindricalChebyshevFftwTransform::generateGrid(const int size, const MHDFloat gapWidth, const MHDFloat rRatio)
   {
      if(gapWidth > 0 && rRatio >= 0)
      {
         // Initialise grid storage
         Array grid(size);

         // Create Chebyshev grid
         for(int k = 0; k < size; k++)
         {
            grid(k) = std::cos((MathConstants::PI)*(static_cast<MHDFloat>(k)+0.5)/static_cast<MHDFloat>(size));

            grid(k) = 0.5*gapWidth*grid(k) + 0.5*gapWidth*(1.+rRatio)/(1-rRatio);
         }

         return grid;
      } else
      {
         throw Exception("generateGrid called with incompatible gap width or radii ratio");
      }
   }

   CylindricalChebyshevFftwTransform::CylindricalChebyshevFftwTransform()
      : mFPlan(NULL), mBPlan(NULL), mGapWidth(-1), mRRatio(-1)
   {
   }

   CylindricalChebyshevFftwTransform::~CylindricalChebyshevFftwTransform()
   {
      // Cleanup memory used by FFTW
      FftwLibrary::cleanupFft();
   }

   void CylindricalChebyshevFftwTransform::init(CylindricalChebyshevFftwTransform::SharedSetupType spSetup)
   {
      // Store the shared pointer to setup object
      this->mspSetup = spSetup;

      // Set the scaling factor
      this->mspSetup->setScale(1.0/static_cast<MHDFloat>(2*this->mspSetup->fwdSize()));

      // Initialise FFTW interface
      this->initFft();

      // Register the FFTW object
      FftwLibrary::registerFft();
   }

   void CylindricalChebyshevFftwTransform::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      list.insert(NonDimensional::GAPWIDTH);

      list.insert(NonDimensional::RRATIO);
   }

   void CylindricalChebyshevFftwTransform::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      this->mGapWidth = options.find(NonDimensional::GAPWIDTH)->second;

      this->mRRatio = options.find(NonDimensional::RRATIO)->second;
   }

   Array CylindricalChebyshevFftwTransform::meshGrid() const
   {
      return CylindricalChebyshevFftwTransform::generateGrid(this->mspSetup->fwdSize(), this->mGapWidth, this->mRRatio);
   }

   void CylindricalChebyshevFftwTransform::initFft()
   {  
      /// \mhdBug implement strideed stranforms for complex <-> complex case if possible

      int fwdSize = this->mspSetup->fwdSize();
      int bwdSize = this->mspSetup->bwdSize();
      int howmany = this->mspSetup->howmany();

      // Create the two plans
      const int  *fftSize = &fwdSize;

      if(this->mspSetup->type() == FftSetup::COMPONENT)
      {
      } else
      {
      }

      // Initialise temporary storage
      this->mTmpIn.setZero(fwdSize, howmany);
      this->mTmpOut.setZero(bwdSize, howmany);

      // Create the physical to spectral plan
      const fftw_r2r_kind fwdKind[] = {FFTW_REDFT10};
      this->mFPlan = fftw_plan_many_r2r(1, fftSize, howmany, this->mTmpIn.data(), NULL, 1, fwdSize, this->mTmpOut.data(), NULL, 1, bwdSize, fwdKind, FftwLibrary::planFlag());

      // Create the spectral to physical plan
      const fftw_r2r_kind bwdKind[] = {FFTW_REDFT01};
      this->mBPlan = fftw_plan_many_r2r(1, fftSize, howmany, this->mTmpOut.data(), NULL, 1, bwdSize, this->mTmpIn.data(), NULL, 1, fwdSize, bwdKind, FftwLibrary::planFlag());
   }

   void CylindricalChebyshevFftwTransform::cleanupFft()
   {
      // Detroy forward plan
      if(this->mFPlan)
      {
         fftw_destroy_plan(this->mFPlan);
      }

      // Detroy backward plan
      if(this->mBPlan)
      {
         fftw_destroy_plan(this->mBPlan);
      }

      // Unregister the FFTW object
      FftwLibrary::unregisterFft();

      // cleanup fftw library
      FftwLibrary::cleanupFft();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   MHDFloat CylindricalChebyshevFftwTransform::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Storage required for the fftw plans 
      mem += 8.0*2.0;

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}
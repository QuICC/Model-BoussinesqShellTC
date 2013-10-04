/** 
 * @file ITauBoundary.cpp
 * @brief Source of the implementation of the spectral boundary conditions
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//

// Class include
//
#include "SpectralOperators/ITauBoundary.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   ITauBoundary::ITauBoundary(const MHDFloat c, const int nN, const int nEq)
      : mC(c), mN(nN), mNeq(nEq), mIsComplex(false), mRTau(0,0), mZTau(0,0)
   {
   }

   ITauBoundary::~ITauBoundary()
   {
   }

   DecoupledZMatrix ITauBoundary::tauLines(const Boundary::BCVector& bcs) const
   {
      assert(this->mNeq >= 0);
      assert(static_cast<unsigned int>(this->mNeq) == bcs.size());

      // Count boundary conditions
      int nBCs = bcs.size();

      // Flags to know if real and imaginary parts are used
      bool hasReal = false;
      bool hasImag = false;

      // Check for compatible sizes
      assert(this->mN >= nBCs);

      // Initialise tau lines matrix
      DecoupledZMatrix lines(this->mN, nBCs);
      lines.setZero();

      // Map iterator
      Boundary::BCVector::const_iterator it;

      // Create boundary values
      int idx = 0;
      for(it = bcs.begin(); it != bcs.end(); ++it)
      {
         switch(it->type)
         {
            case Boundary::VALUE:
               lines.real().col(idx) = this->value(it->position);
               idx++;
               hasReal = true;
               break;
            case Boundary::D1:
               lines.real().col(idx) = this->firstDerivative(it->position);
               idx++;
               hasReal = true;
               break;
            case Boundary::D2:
               lines.real().col(idx) = this->secondDerivative(it->position);
               idx++;
               hasReal = true;
               break;
            case Boundary::BETA_SLOPE:
               /// \warning Beta slope boundary conditions does NOT include the factors
               if(it->position == Boundary::RIGHT)
               {
                  lines.imag().col(idx) = static_cast<MHDFloat>(-1)*this->value(it->position);
               } else //if(it->position == Boundary::LEFT)
               {
                  lines.imag().col(idx) = this->value(it->position);
               }
               idx++;
               hasImag = true;
               break;
            default:
               throw Exception("Unknown boundary condition ID for tau method");
               break;
         }
      }

      // Clear real matrix if not conditions have been set
      if(!hasReal)
      {
         lines.real().resize(0,0);
      }

      // Clear imaginary matrix if not conditions have been set
      if(!hasImag)
      {
         lines.imag().resize(0,0);
      }

      return lines;
   }

   void ITauBoundary::createTauMatrix(const Boundary::BCVector& bcs)
   {
      // Builde the tau lines
      DecoupledZMatrix lines = ITauBoundary::tauLines(bcs);

      // Create sparse matrices
      if(lines.imag().size() != 0)
      {
         this->mZTau.resize(this->mN, this->mN);
         this->mZTau.reserve(lines.real().rows()*lines.real().cols() + lines.imag().rows()*lines.imag().cols());
         this->mIsComplex = true;
      } else
      {
         this->mRTau.resize(this->mN, this->mN);
         this->mRTau.reserve(lines.real().rows()*lines.real().cols());
         this->mIsComplex = false;
      }

      // Fill complex Tau line matrix
      if(this->mIsComplex)
      {
         // Make complex matrix
         MatrixZ tmp(this->mZTau.cols(), lines.imag().cols());
         tmp.setZero();
         if(lines.real().size() != 0)
         {
            tmp.real() = lines.real();
         }
         tmp.imag() = lines.imag();

         for(int j = 0; j < this->mZTau.cols(); ++j)
         {
            // create column j
            this->mZTau.startVec(j);

            for(int i = 0; i < tmp.cols(); ++i)
            {
               this->mZTau.insertBack(i,j) = this->mC*tmp(j,i);
            }
         }
         this->mZTau.prune(MHDComplex(1e-32,1e-32));
         this->mZTau.finalize(); 

      // Fill real Tau line matrix
      } else
      {
         for(int j = 0; j < this->mRTau.cols(); ++j)
         {
            // create column j
            this->mRTau.startVec(j);

            for(int i = 0; i < lines.real().cols(); ++i)
            {
               this->mRTau.insertBack(i,j) = this->mC*lines.real()(j,i);
            }
         }
         this->mRTau.prune(1e-32);
         this->mRTau.finalize(); 
      }
   }

}
}

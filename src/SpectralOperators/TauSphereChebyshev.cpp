/** 
 * @file TauSphereChebyshev.cpp
 * @brief Source of the implementation of the spectral Chebyshev Tau operators for a sphere radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//
#include <cassert>

// External includes
//
#include <Eigen/Sparse>

// Class include
//
#include "SpectralOperators/TauSphereChebyshev.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   TauSphereChebyshev::TauSphereChebyshev(const MHDFloat c, const int nN, const Boundary::BCVector& bcs, const int nEq)
      : ITauBoundary(c, nN, nEq)
   {
      this->createTauMatrix(bcs);
   }

   TauSphereChebyshev::~TauSphereChebyshev()
   {
   }

   Array TauSphereChebyshev::value(Boundary::BCPosition pt) const
   {
      if(pt == Boundary::LEFT)
      {
         return this->leftUnit();
      } else
      {
         return this->rightUnit();
      }
   }

   Array TauSphereChebyshev::firstDerivative(Boundary::BCPosition pt) const
   {
      Array tau = Array::Zero(this->mN);
      for(int i = 1; i < tau.size(); i++)
      {
         tau(i) = static_cast<MHDFloat>(i*i);
      }

      if(pt == Boundary::LEFT)
      {
         return tau.array()*(-1.0*this->leftUnit().array());
      }else 
      {
         return tau.array()*this->rightUnit().array();
      }
   }

   Array TauSphereChebyshev::secondDerivative(Boundary::BCPosition pt) const
   {
      Array tau = Array::Zero(this->mN);

      for(int i = 2; i < tau.size(); i++)
      {
         tau(i) = (1.0/3.0)*static_cast<MHDFloat>(std::pow(i,4) - std::pow(i,2));
      }

      if(pt == Boundary::LEFT)
      {
         return tau.array()*this->leftUnit().array();
      } else
      {
         return tau.array()*this->rightUnit().array();
      }
   }

   Array TauSphereChebyshev::leftUnit() const
   {
      Array unit(this->mN);
      for(int i = 0; i < unit.size(); i++)
      {
         unit(i) = this->c(i)*std::pow(-1.0,i);
      }

      return unit;
   }

   Array TauSphereChebyshev::rightUnit() const
   {
      Array unit(this->mN);
      for(int i = 0; i < unit.size(); i++)
      {
         unit(i) = this->c(i);
      }

      return unit;
   }

   MHDFloat TauSphereChebyshev::c(const int n) const
   {
      if(n == 0)
      {
         return 1.0;

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         #ifdef GEOMHDISCC_CHEBYSHEV_HAS_C
            return 1.0;
         #else
            return 2.0;
         #endif
      }
   }

   MHDFloat TauSphereChebyshev::c_1(const int n) const
   {
      if(n == 0)
      {
         return 1.0;

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         #ifdef GEOMHDISCC_CHEBYSHEV_HAS_C
            return 1.0;
         #else
            return 0.5;
         #endif
      }
   }

}
}
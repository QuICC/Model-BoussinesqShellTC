/** \file ChebyshevBoundary.cpp
 *  \brief Source of the implementation of the spectral boundary operator for the Chebyshev basis
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "SpectralOperators/ChebyshevBoundary.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   ChebyshevBoundary::ChebyshevBoundary(const int basisN)
      : IBoundary(basisN)
   {
   }

   ChebyshevBoundary::~ChebyshevBoundary()
   {
   }

   MHDFloat ChebyshevBoundary::c(const int n) const
   {
      if(n == 0)
      {
         return 2.0;

      } else if(n < 0)
      {
         return 0.0;

      } else
      {
         return 1.0;
      }
   }

   Array ChebyshevBoundary::leftUnit() const
   {
      Array unit(this->basisN());
      for(int i = 0; i < unit.size(); i++)
      {
         unit(i) = (1.0/this->c(i))*std::pow(-1.0,i);
      }
   }

   Array ChebyshevBoundary::rightUnit() const
   {
      Array unit(this->basisN());
      for(int i = 0; i < unit.size(); i++)
      {
         unit(i) = (1.0/this->c(i));
      }
   }

   Array ChebyshevBoundary::value(IBoundary::Position pt) const
   {
      if(pt == LEFT)
      {
         return this->leftUnit();
      } else
      {
         return this->rightUnit();
      }
   }

   Array ChebyshevBoundary::firstDerivative(IBoundary::Position pt) const
   {
      Array tau = Array::zero(this->basisN());
      for(int i = 1; i < val.size(); i++)
      {
         tau(i) = static_cast<MHDFloat>(i*i);
      }

      if(pt == LEFT)
      {
         return tau.array()*(-1.0*this->leftUnit().array());
      }else 
      {
         return tau.array()*this->rightUnit().array();
      }
   }

   Array ChebyshevBoundary::secondDerivative(IBoundary::Position pt) const
   {
      Array tau = Array::zero(this->basisN());

      for(int i = 2; i < val.size(); i++)
      {
         tau(i) = (1.0/3.0)*static_cast<MHDFloat>(std::pow(i,4) - std::pow(i,2));
      }

      if(pt == LEFT)
      {
         return tau.array()*this->leftUnit().array();
      } else
      {
         return tau.array()*this->rightUnit().array();
      }
   }

}
}
/** \file ScalarFieldSetup.cpp
 *  \brief Source of the setup object for the scalar fields
 */

// System includes
//

// External includes
//

// Class include
//
#include "ScalarFields/ScalarFieldSetup.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Datatypes {

   ScalarFieldSetup::ScalarFieldSetup(const std::vector<ArrayI>& idx1D, const std::vector<ArrayI>& idx2D, const ArrayI& idx3D)
   {
      // Get first dimension size (not fully general case as sliced can't deal with it)
      this->mspDim1D = SharedArrayI(new ArrayI(idx1D.size()));
      for(int j=0;  j < idx1D.size(); j++)
      {
         (*this->mspDim1D)(j) = idx1D.at(j).rows();
      }

      // Get second dimension size
      this->mspDim2D = SharedArrayI(new ArrayI(idx2D.size()));
      for(int j=0;  j < idx2D.size(); j++)
      {
         (*this->mspDim2D)(j) = idx2D.at(j).size();
      }

      // Get third dimension size
      this->mDim3D = idx3D.size();
   }

   ScalarFieldSetup::ScalarFieldSetup(SharedArrayI spDim1D, SharedArrayI spDim2D, const int dim3D)
      : mspDim1D(spDim1D), mspDim2D(spDim2D), mDim3D(dim3D)
   {
   }

   ScalarFieldSetup::~ScalarFieldSetup()
   {
   }

}
}

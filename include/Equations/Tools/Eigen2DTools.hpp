/**
 * @file Eigen2DTools.hpp
 * @brief Implementation of some tools for schemes with two eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EIGEN2DTOOLS_HPP
#define EIGEN2DTOOLS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include<vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Equations/Tools/IEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Tools for equations with two eigen directions
    */
   class Eigen2DTools: public IEigenTools
   {
      public:
         /**
          * @brief Constructor
          */
         Eigen2DTools();

         /**
          * @brief Destructor
          */
         virtual ~Eigen2DTools();

      private:

         /**
          * @brief Set eigen values
          */
         virtual std::vector<MHDFloat> identifyEigs(const SharedResolution& spRes, const int matIdx) const;

         /**
          * @brief Compute the number of matrix operators for field coupling
          *
          * @param spRes   Shared resolution
          */
         virtual int computeNMat(const SharedResolution& spRes) const;

         /**
          * @brief Interpret Tau resolution provided by python code
          */
         virtual void interpretTauN(ArrayI& rTauNs, const SharedResolution& spRes) const;

         /**
          * @brief Interpret Galerkin resolution provided by python code
          */
         virtual void interpretGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const;

         /**
          * @brief Interpret number of RHS provided by python code
          */
         virtual void interpretRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const;

         /**
          * @brief Interpret system size provided by python code
          */
         virtual void interpretSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const;
   };

   /// Typedef for a shared Eigen2DTools
   typedef SharedPtrMacro<Eigen2DTools> SharedEigen2DTools;

}
}

#endif // EIGEN2DTOOLS_HPP

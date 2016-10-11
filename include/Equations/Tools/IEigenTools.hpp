/**
 * @file IEigenTools.hpp
 * @brief Interface for the equations eigen direction tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IEIGENTOOLS_HPP
#define IEIGENTOOLS_HPP

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
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Interface for equations eigen directions tools
    */
   class IEigenTools
   {
      public:
         /**
          * @brief Constructor
          */
         IEigenTools();

         /**
          * @brief Destructor
          */
         virtual ~IEigenTools();

         /**
          * @brief Get eigen values
          */
         std::vector<MHDFloat> getEigs(const SharedResolution& spRes, const int matIdx) const;

         /**
          * @brief Get the number of matrix operators for field coupling
          *
          * @param spRes   Shared resolution
          */
         int nMat(const SharedResolution& spRes) const;

         /**
          * @brief Set Tau resolution
          */
         void setTauN(ArrayI& rTauNs, const SharedResolution& spRes) const;

         /**
          * @brief Set Galerkin resolution provided by python code
          */
         void setGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const;

         /**
          * @brief Set number of RHS provided by python code
          */
         void setRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const;

         /**
          * @brief Set system size provided by python code
          */
         void setSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const;

      private:
         /**
          * @brief Create set of eigen values
          */
         virtual std::vector<MHDFloat> identifyEigs(const SharedResolution& spRes, const int matIdx) const = 0;

         /**
          * @brief Compute the number of matrix operators for field coupling
          *
          * @param spRes   Shared resolution
          */
         virtual int computeNMat(const SharedResolution& spRes) const = 0;

         /**
          * @brief Interpret Tau resolution provided by python code
          */
         virtual void interpretTauN(ArrayI& rTauNs, const SharedResolution& spRes) const = 0;

         /**
          * @brief Interpret Galerkin resolution provided by python code
          */
         virtual void interpretGalerkinN(ArrayI& rGalerkinNs, const SharedResolution& spRes) const = 0;

         /**
          * @brief Interpret number of RHS provided by python code
          */
         virtual void interpretRhsN(ArrayI& rRhsCols, const SharedResolution& spRes) const = 0;

         /**
          * @brief Interpret system size provided by python code
          */
         virtual void interpretSystemN(ArrayI& rSystemNs, const SharedResolution& spRes) const = 0;
   };

   /// Typedef for a shared IEigenTools
   typedef SharedPtrMacro<IEigenTools> SharedIEigenTools;

}
}

#endif // IEIGENTOOLS_HPP

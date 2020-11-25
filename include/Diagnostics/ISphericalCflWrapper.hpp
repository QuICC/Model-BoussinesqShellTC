/**
 * @file ISphericalCflWrapper.hpp
 * @brief CFL constraint in a spherical geometry 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP
#define QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/NonDimensional.hpp"
#include "Diagnostics/ICflWrapper.hpp"
#include "Diagnostics/IVectorWrapper.hpp"

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical geometry
    */
   class ISphericalCflWrapper: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          */
         ISphericalCflWrapper(const SharedIVectorWrapper spVelocity, const std::map<NonDimensional::Id,MHDFloat>& params);

         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          */
         ISphericalCflWrapper(const SharedIVectorWrapper spVelocity, const SharedIVectorWrapper spMagnetic, const std::map<NonDimensional::Id,MHDFloat>& params);

         /**
          * @brief Constructor
          */
         ~ISphericalCflWrapper();

         /**
          * @brief Initialize wrapper
          */
         virtual void init(const std::vector<Array>& mesh);

         /**
          * @brief Get initial CFL constraint
          */
         virtual Matrix initialCfl() const;

         /**
          * @brief Get CFL constraint
          */
         virtual Matrix cfl() const;

      protected:

      private:
         /**
          * @brief Get effective max harmonic degree L
          */
         virtual MHDFloat effectiveMaxL(const MHDFloat r) const = 0;

         /**
          * @brief Initialise the mesh spacings
          */
         void initMesh(const std::vector<Array>& mesh);

         /**
          * @brief Courant constant used for the CFL computation
          */
         const MHDFloat mcCourant;

         /**
          * @brief Alfven wave scale
          */
         const MHDFloat mcAlfvenScale;

         /**
          * @brief Alfven wave damping
          */
         const MHDFloat mcAlfvenDamping;

         /**
          * @brief CFL conditions
          */
         Array mGlobalCfl;

         /**
          * @brief Spacing between grid points
          */
         std::vector<Array> mMeshSpacings;
   };

   /// Typedef for a shared ISphericalCflWrapper
   typedef SharedPtrMacro<ISphericalCflWrapper> SharedISphericalCflWrapper;
}
}

#endif // QUICC_DIAGNOSTICS_ISPHERICALCFLWRAPPER_HPP

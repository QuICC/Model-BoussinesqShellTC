/**
 * @file SphericalCflWrapper.hpp
 * @brief CFL constraint in a spherical geometry 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERICALCFLWRAPPER_HPP
#define SPHERICALCFLWRAPPER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Diagnostics/ICflWrapper.hpp"
#include "Diagnostics/IVelocityWrapper.hpp"

namespace GeoMHDiSCC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a spherical geometry
    */
   class SphericalCflWrapper: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          * @param Vector of physical space grid
          */
         SphericalCflWrapper(const SharedIVelocityWrapper spVelocity, const std::vector<Array>& mesh);

         /**
          * @brief Constructor
          */
         ~SphericalCflWrapper();

         /**
          * @brief Get initial CFL constraint
          */
         virtual MHDFloat initialCfl() const;

         /**
          * @brief Get CFL constraint
          */
         virtual MHDFloat cfl() const;

      protected:

      private:
         /**
          * @brief Initialise the mesh spacings
          */
         void initMesh(const std::vector<Array>& mesh);

         /**
          * @brief Courant constant used for the CFL computation
          */
         const MHDFloat mcCourant;

         /**
          * @brief Spacing between grid points
          */
         std::vector<Array> mMeshSpacings;
   };

   /// Typedef for a shared SphericalCflWrapper
   typedef SharedPtrMacro<SphericalCflWrapper> SharedSphericalCflWrapper;
}
}

#endif // SPHERICALCFLWRAPPER_HPP

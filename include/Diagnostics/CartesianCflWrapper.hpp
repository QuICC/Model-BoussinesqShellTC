/**
 * @file CartesianCflWrapper.hpp
 * @brief CFL constraint in a Cartesian geometry 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIANCFLWRAPPER_HPP
#define CARTESIANCFLWRAPPER_HPP

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

namespace QuICC {

namespace Diagnostics {

   /**
    * @brief CFL constraint in a Cartesian geometry
    */
   class CartesianCflWrapper: public ICflWrapper
   {
      public:
         /**
          * @brief Constructor
          *
          * @param Velocity wrapper
          */
         CartesianCflWrapper(const SharedIVelocityWrapper spVelocity);

         /**
          * @brief Constructor
          */
         ~CartesianCflWrapper();

         /**
          * @brief Initialize wrapper
          */
         virtual void init(const std::vector<Array>& mesh);

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

   /// Typedef for a shared CartesianCflWrapper
   typedef SharedPtrMacro<CartesianCflWrapper> SharedCartesianCflWrapper;
}
}

#endif // CARTESIANCFLWRAPPER_HPP

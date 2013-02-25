/** \file ScalarFieldSetup.hpp
 *  \brief Single configuration class for the different scalar fields
 */

#ifndef SCALARFIELDSETUP_HPP
#define SCALARFIELDSETUP_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Single configuration class for the different scalar fields
    */
   class ScalarFieldSetup
   {
      public:
         /**
          * @brief Constructor based on index information
          *
          * This constructor acutally allocates the memory for the shared data
          *
          * @param idx1D  Index information of the first dimension 
          * @param idx2D  Index information of the second dimension 
          * @param idx3D  Index information of the third dimension 
          */
         ScalarFieldSetup(const std::vector<ArrayI>& idx1D, const std::vector<ArrayI>& idx2D, const ArrayI& idx3D);

         /**
          * @brief Constructor
          *
          * @param spDim1D  Size of the first dimension 
          * @param spDim2D  Size of the second dimension 
          * @param dim3D  Size of the third dimension 
          */
         ScalarFieldSetup(SharedArrayI spDim1D, SharedArrayI spDim2D = SharedArrayI(), const int dim3D = 0);

         /**
          * @brief Destructor
          */
         virtual ~ScalarFieldSetup();

         /**
          * @brief Get first dimension
          */
         SharedArrayI spDim1D() const;

         /**
          * @brief Get second dimension
          */
         SharedArrayI spDim2D() const;

         /**
          * @brief Get third dimension
          */
         int dim3D() const;
         
      protected:
         /**
          * @brief Size in the first dimension
          */
         SharedArrayI mspDim1D;

         /**
          * @brief Size in the second dimension
          */
         SharedArrayI mspDim2D;

         /**
          * @brief Size in the third dimension
          */
         int mDim3D;

      private:
   };

   inline SharedArrayI ScalarFieldSetup::spDim1D() const
   {
      return this->mspDim1D;
   }

   inline SharedArrayI ScalarFieldSetup::spDim2D() const
   {
      return this->mspDim2D;
   }

   inline int ScalarFieldSetup::dim3D() const
   {
      return this->mDim3D;
   }

   /// Typedef for an smart reference counting pointer for a ScalarFieldSetup
   typedef SharedPtrMacro<ScalarFieldSetup>   SharedScalarFieldSetup;

}
}

#endif // SCALARFIELDSETUP_HPP

/** \file VectorField.hpp
 *  \brief Implementation of a generic vector field
 *
 *  \mhdBug Needs test
 */

#ifndef VECTORFIELD_HPP
#define VECTORFIELD_HPP

// Configuration includes
//

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
    * @brief Implementation of a generic vector field
    */
   template <typename TScalar, int COMPONENTS> class VectorField
   {
      public:
         /**
          * @brief Constructor with simplified interface
          *
          * @param setup Setup object for the scalar fields
          */
         VectorField(const typename TScalar::SetupType& setup);

         /**
          * @brief Destructor
          */
         virtual ~VectorField();

         /**
          * @brief Get field component
          *
          * @param i Index of the component
          */
         const TScalar& comp(const int i) const;

         /**
          * @brief Set field component
          *
          * @param i Index of the component
          */
         TScalar& rComp(const int i);

         /**
          * @brief Get field components
          */
         const std::vector<TScalar>& data() const;

         /**
          * @brief Set field components
          */
         std::vector<TScalar>& rData();

         /**
          * @brief Initialise field components to zero
          */
         void initialiseZeros();

      #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:
         /**
          * @brief Storage for the field components
          */
         std::vector<TScalar> mComponents;

      private:
   };

   template <typename TScalar, int COMPONENTS> inline const TScalar& VectorField<TScalar,COMPONENTS>::comp(const int i) const
   {
      // Assert that index is valid
      assert(this->mComponents.size() > i);

      return this->mComponents.at(i);
   }

   template <typename TScalar, int COMPONENTS> inline TScalar& VectorField<TScalar,COMPONENTS>::rComp(const int i)
   {
      // Assert that index is valid
      assert(this->mComponents.size() > i);

      return this->mComponents.at(i);
   }

   template <typename TScalar, int COMPONENTS> VectorField<TScalar, COMPONENTS>::VectorField(const typename TScalar::SetupType& setup)
   {
      // Initialise the components
      for(int i = 0; i < COMPONENTS; i++)
      {
         this->mComponents.push_back(TScalar(setup));
      }
   }

   template <typename TScalar, int COMPONENTS> VectorField<TScalar, COMPONENTS>::~VectorField()
   {
   }

   template <typename TScalar, int COMPONENTS> inline const std::vector<TScalar>& VectorField<TScalar,COMPONENTS>::data() const
   {
      return this->mComponents;
   }

   template <typename TScalar, int COMPONENTS> inline std::vector<TScalar>& VectorField<TScalar,COMPONENTS>::rData()
   {
      return this->mComponents;
   }

   template <typename TScalar, int COMPONENTS> void VectorField<TScalar, COMPONENTS>::initialiseZeros()
   {
      // Initialise the components
      typename std::vector<TScalar>::iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         it->initialiseZeros();
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TScalar, int COMPONENTS> MHDFloat VectorField<TScalar, COMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;
      typename std::vector<TScalar>::const_iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         mem += it->requiredStorage();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // VECTORFIELD_HPP

/** \file VariableDomain.hpp
 *  \brief Implementation of the variable's domain
 *
 *  \mhdBug Needs test
 */

#ifndef VARIABLEDOMAIN_HPP
#define VARIABLEDOMAIN_HPP

// Configuration includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * \brief Implementation for a full sphere field
    *
    * \tparam TVariable Type of the variable
    * \tparam DOMAINS   Number of different domains
    */
   template <typename TVariable, int DOMAINS> class VariableDomain
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spRes Resolution information
          */
         VariableDomain(SharedResolution spRes);

         /**
          * @brief Destructor
          */
         virtual ~VariableDomain();

         /**
          * @brief Get variable of requested domain
          *
          * @param i Index of the domain
          */
         const TVariable&  dom(const int i) const;

         /**
          * @brief Set Physical variable in full sphere
          *
          * @param i Index of the domain
          */
         TVariable&  rDom(const int i);

         /**
          * @brief initialise to zeros
          */
         void setZeros();

         /**
          * @brief Initialise the physical values storage
          */
         void initPhysical();

         /**
          * @brief Initialise the physical differential values storage
          */
         void initPhysicalDiff();

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:
         /**
          * @brief Storage for the OC (or full sphere) part
          */
         std::vector<TVariable>   mDomains;

      private:
   };

   template <typename TVariable, int DOMAINS> inline const TVariable& VariableDomain<TVariable,DOMAINS>::dom(const int i) const
   {
      // Assert for valid index
      assert(this->mDomains.size() > static_cast<size_t>(i));

      return this->mDomains.at(i);
   }

   template <typename TVariable, int DOMAINS> inline TVariable& VariableDomain<TVariable,DOMAINS>::rDom(const int i)
   {
      // Assert for valid index
      assert(i >= 0);
      assert(this->mDomains.size() > static_cast<size_t>(i));

      return this->mDomains.at(i);
   }

   template <typename TVariable, int DOMAINS> VariableDomain<TVariable,DOMAINS>::VariableDomain(SharedResolution spRes)
   {
      for(int i = 0; i < DOMAINS; i++)
      {
         this->mDomains.push_back(TVariable(spRes));
      }
   }

   template <typename TVariable, int DOMAINS> VariableDomain<TVariable,DOMAINS>::~VariableDomain()
   {
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::setZeros()
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).setZeros();
      }
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::initPhysical()
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).initPhysical();
      }
   }

   template <typename TVariable, int DOMAINS> void  VariableDomain<TVariable,DOMAINS>::initPhysicalDiff()
   {
      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         this->mDomains.at(i).initPhysicalDiff();
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TVariable, int DOMAINS> MHDFloat  VariableDomain<TVariable,DOMAINS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Loop over all domains
      for(size_t i = 0; i < this->mDomains.size(); i++)
      {
         mem += this->mDomains.at(i).requiredStorage();
      }
      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // VARIABLEDOMAIN_HPP

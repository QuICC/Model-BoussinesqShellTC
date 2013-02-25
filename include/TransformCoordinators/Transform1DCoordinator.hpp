/** \file Transform1DCoordinator.hpp
 *  \brief Implementation of a transform coordinator in 1D
 */

#ifndef TRANSFORM1DCOORDINATOR_HPP
#define TRANSFORM1DCOORDINATOR_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/PhysicalNames.hpp"
#include "Resolutions/Resolution.hpp"


namespace GeoMHDiSCC {

   /**
    * 
    * @brief Implementation of a transform coordinator in 1D
    *
    * The transform coordinator overlooks the whole transform process. Providing storage, communication
    * and the actual transforms
    */ 
   template <typename T1D, typename TCommunicator> class Transform1DCoordinator
   {
      public:
         /// Typedef for the first dimension transform
         typedef T1D Transform1DType;

         /// Typedef for the communicator type
         typedef TCommunicator CommunicatorType;

         /**
          * @brief Very basic constructor
          */
         Transform1DCoordinator();

         /**
          * @brief Destructor
          */
         virtual ~Transform1DCoordinator();

         /**
          * @brief Initialise transform of the 1D coordinator
          *
          * @param setup1D Setup object for 1D transform
          * @param varInfo Variables information
          */
         void initTransforms(const typename T1D::SetupType& setup1D, const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Initialise the data communicator
          *
          * @param spRes   Resolution information object
          */
         void initCommunicator(SharedResolution spRes);

         /**
          * @brief Get the transform for the first dimension
          */
         T1D&   transform1D();

         /**
          * @brief Get the communicator
          */
         TCommunicator&  communicator();

         /**
          * @brief Need to compute backward transform to physical values?
          *
          * @param name Name of the field
          */
         bool needPhysical(PhysicalNames::Id name);

         /**
          * @brief Need to compute backward transform to physical differential values (gradient or curl)?
          */
         bool needPhysicalDiff(PhysicalNames::Id name);

         /**
          * @brief Get grid array(s) of the mesh
          */
         std::vector<Array>  mesh();
         
      protected:

         /**
          * @brief The communicator
          */
         TCommunicator  mCommunicator;

      private:
         /**
          * @brief The first dimension transform
          */
         T1D  mTransform1D;

         /**
          * @brief Variable information
          */
         std::map<PhysicalNames::Id, TriBool>  mVarInfo;

         /**
          * @brief Initialise the transforms
          *
          * @param setup1D Setup object for 1D transform
          */
         void initTransform(const typename T1D::SetupType& setup1D);
   };

   template <typename T1D, typename TCommunicator> inline T1D&  Transform1DCoordinator<T1D,TCommunicator>::transform1D()
   {
      return this->mTransform1D;
   }

   template <typename T1D, typename TCommunicator> inline TCommunicator&  Transform1DCoordinator<T1D,TCommunicator>::communicator()
   {
      return this->mCommunicator;
   }

   template <typename T1D, typename TCommunicator> inline bool Transform1DCoordinator<T1D,TCommunicator>::needPhysical(PhysicalNames::Id name)
   {
      // Safety assertion
      assert(this->mVarInfo.find(name) != this->mVarInfo.end());

      return this->mVarInfo.find(name)->second(1);
   }

   template <typename T1D, typename TCommunicator> inline bool Transform1DCoordinator<T1D,TCommunicator>::needPhysicalDiff(PhysicalNames::Id name)
   {
      // Safety assertion
      assert(this->mVarInfo.find(name) != this->mVarInfo.end());

      return this->mVarInfo.find(name)->second(2);
   }

   template <typename T1D, typename TCommunicator> Transform1DCoordinator<T1D,TCommunicator>::Transform1DCoordinator()
   {
   }

   template <typename T1D, typename TCommunicator> Transform1DCoordinator<T1D,TCommunicator>::~Transform1DCoordinator()
   {
   }

   template <typename T1D, typename TCommunicator> void Transform1DCoordinator<T1D,TCommunicator>::initTransforms(const typename T1D::SetupType& setup1D, const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {
      // Initialise the transforms
      this->initTransform(setup1D);

      // Store information about variables
      std::map<PhysicalNames::Id, std::pair<bool,TriBool> >::const_iterator   it;
      for(it = varInfo.begin(); it != varInfo.end(); it++)
      {
         this->mVarInfo.insert(std::make_pair(it->first, it->second.second));
      }
   }

   template <typename T1D, typename TCommunicator> void Transform1DCoordinator<T1D,TCommunicator>::initTransform(const typename T1D::SetupType& setup1D)
   {
      // initialise the first transform
      this->mTransform1D.init(setup1D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem1D = this->mTransform1D.requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORMS, mem1D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM1D, mem1D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

  template <typename T1D, typename TCommunicator> void Transform1DCoordinator<T1D,TCommunicator>::initCommunicator(SharedResolution spRes)
   {
      // initialise the communicator
      this->mCommunicator.init(*spRes->spFwdSetup(0), *spRes->spBwdSetup(0));
   }

   template <typename T1D, typename TCommunicator> std::vector<Array> Transform1DCoordinator<T1D,TCommunicator>::mesh()
   {
      std::vector<Array>   mesh;

      // Add grid of first dimension
      mesh.push_back(this->transform1D().meshGrid());

      return mesh;
   }

}

#endif // TRANSFORM1DCOORDINATOR_HPP

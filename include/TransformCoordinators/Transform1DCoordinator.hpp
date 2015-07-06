/** 
 * @file Transform1DCoordinator.hpp
 * @brief Implementation of a transform coordinator in 1D
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORM1DCOORDINATOR_HPP
#define TRANSFORM1DCOORDINATOR_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//
#include "TypeSelectors/TransformTreeSelector.hpp"

// System includes
//
#include<set>
#include<map>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Enums/NonDimensional.hpp"
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
          * @param spRes            Resolution information object
          * @param integratorTree   Transform integrator tree
          * @param projectorTree    Transform projector tree
          */
         void initTransforms(SharedResolution spRes, const std::vector<Transform::IntegratorTree>& integratorTree, const std::vector<Transform::ProjectorTree>& projectorTree);

         /**
          * @brief Initialise the data communicator
          *
          * @param spRes   Resolution information object
          */
         void initCommunicator(SharedResolution spRes);

         /**
          * @brief Get list of required options for transform
          */
         void requiredOptions(std::set<NonDimensional::Id>& list) const;

         /**
          * @brief Set options of transform
          */
         void setOptions(const std::map<NonDimensional::Id, MHDFloat>& options);

         /**
          * @brief Get the transform for the first dimension
          */
         T1D&   transform1D();

         /**
          * @brief Get the communicator
          */
         TCommunicator&  communicator();

         /**
          * @brief Get the transform integrator tree
          */
         const std::vector<Transform::IntegratorTree>& integratorTree() const;

         /**
          * @brief Get the transform projector tree
          */
         const std::vector<Transform::ProjectorTree>& projectorTree() const;

         /**
          * @brief Need to compute backward transform to physical values?
          *
          * @param name Name of the field
          */
         bool needPhysical(PhysicalNames::Id name);

         /**
          * @brief Need to compute backward transform to physical gradient values (scalar or vector)?
          */
         bool needPhysicalGradient(PhysicalNames::Id name);

         /**
          * @brief Need to compute backward transform to physical curl values?
          */
         bool needPhysicalCurl(PhysicalNames::Id name);

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
          * @brief Transform integrator tree
          */
         std::vector<Transform::IntegratorTree> mIntegratorTree;

         /**
          * @brief Transform projector tree
          */
         std::vector<Transform::ProjectorTree> mProjectorTree;

         /**
          * @brief Initialise the transforms
          *
          * @param spSetup1D Shared setup object for 1D transform
          */
         void initTransform(typename T1D::SharedSetupType spSetup1D);
   };

   template <typename T1D, typename TCommunicator> inline T1D&  Transform1DCoordinator<T1D,TCommunicator>::transform1D()
   {
      return this->mTransform1D;
   }

   template <typename T1D, typename TCommunicator> inline TCommunicator&  Transform1DCoordinator<T1D,TCommunicator>::communicator()
   {
      return this->mCommunicator;
   }

   template <typename T1D, typename TCommunicator> inline const std::vector<Transform::IntegratorTree>& Transform1DCoordinator<T1D,TCommunicator>::integratorTree() const
   {
      return this->mIntegratorTree;
   }

   template <typename T1D, typename TCommunicator> inline const std::vector<Transform::ProjectorTree>& Transform1DCoordinator<T1D,TCommunicator>::projectorTree() const
   {
      return this->mProjectorTree;
   }

   template <typename T1D, typename TCommunicator> Transform1DCoordinator<T1D,TCommunicator>::Transform1DCoordinator()
   {
   }

   template <typename T1D, typename TCommunicator> Transform1DCoordinator<T1D,TCommunicator>::~Transform1DCoordinator()
   {
   }

   template <typename T1D, typename TCommunicator> void Transform1DCoordinator<T1D,TCommunicator>::initTransforms(SharedResolution spRes, const std::vector<Transform::IntegratorTree>& integratorTree, const std::vector<Transform::ProjectorTree>& projectorTree)
   {
      // Initialise the transforms
      this->initTransform(std::tr1::static_pointer_cast<typename T1D::SetupType>(spRes->spTransformSetup(Dimensions::Transform::TRA1D)));

      // Store the projector tree
      this->mIntegratorTree = integratorTree;

      // Store the projector tree
      this->mProjectorTree = projectorTree;
   }

   template <typename T1D, typename TCommunicator> void Transform1DCoordinator<T1D,TCommunicator>::initTransform(typename T1D::SharedSetupType spSetup1D)
   {
      // initialise the first transform
      this->mTransform1D.init(spSetup1D);

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
      this->mCommunicator.init(spRes->spFwdSetup(Dimensions::Transform::TRA1D), spRes->spBwdSetup(Dimensions::Transform::TRA1D));
   }

  template <typename T1D, typename TCommunicator> void Transform1DCoordinator<T1D,TCommunicator>::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      this->mTransform1D.requiredOptions(list, Dimensions::Transform::TRA1D);
   }

  template <typename T1D, typename TCommunicator> void Transform1DCoordinator<T1D,TCommunicator>::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      this->mTransform1D.setOptions(options, Dimensions::Transform::TRA1D);
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

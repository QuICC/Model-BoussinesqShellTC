/** 
 * @file Transform3DCoordinator.hpp
 * @brief Implementation of a transform coordinator in 3D
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORM3DCOORDINATOR_HPP
#define TRANSFORM3DCOORDINATOR_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//
#include<set>
#include<map>

// External includes
//

// Project includes
//
#include "Enums/NonDimensional.hpp"
#include "TransformCoordinators/Transform2DCoordinator.hpp"

namespace GeoMHDiSCC {

   /**
    * 
    * @brief Implementation of a transform coordinator in 3D
    *
    * The transform coordinator overlooks the whole transform process. Providing storage, communication
    * and the actual transforms
    */ 
   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> class Transform3DCoordinator: public Transform2DCoordinator<T1D, T2D, TCommunicator>
   {
      public:
         /// Typedef for the first dimension transform
         typedef T1D Transform1DType;

         /// Typedef for the second dimension transform
         typedef T2D Transform2DType;

         /// Typedef for the third dimension transform
         typedef T3D Transform3DType;

         /// Typedef for the communicator type
         typedef TCommunicator CommunicatorType;

         /**
          * @brief Very basic constructor
          */
         Transform3DCoordinator();

         /**
          * @brief Destructor
          */
         virtual ~Transform3DCoordinator();

         /**
          * @brief Initialise the transforms of the 3D coordinator
          *
          * @param spRes   Resolution information object
          * @param varInfo Variables information
          * @param projectorTree Transform projector tree
          */
         void initTransforms(SharedResolution spRes, const std::vector<Transform::ProjectorTree>& projectorTree);

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
         T3D&   transform3D();

         /**
          * @brief Get grid array(s) of the mesh
          */
         std::vector<Array>  mesh();
         
      protected:

      private:
         /**
          * @brief The third dimension transform
          */
         T3D  mTransform3D;

         /**
          * @brief Initialise the transforms
          *
          * @param spSetup3D Shared setup object for 3D transform
          */
         void initTransform(typename T3D::SharedSetupType spSetup3D);
   };

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> inline T3D&  Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::transform3D()
   {
      return this->mTransform3D;
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::Transform3DCoordinator()
      : Transform2DCoordinator<T1D,T2D, TCommunicator>()
   {
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::~Transform3DCoordinator()
   {
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::initTransforms(SharedResolution spRes, const std::vector<Transform::ProjectorTree>& projectorTree)
   {
      // initialise 2 other dimensions
      Transform2DCoordinator<T1D,T2D, TCommunicator>::initTransforms(spRes, projectorTree);

      // Initialise the transforms
      this->initTransform(std::tr1::static_pointer_cast<typename T3D::SetupType>(spRes->spTransformSetup(Dimensions::Transform::TRA3D)));
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::initTransform(typename T3D::SharedSetupType spSetup3D)
   {
      // initialise the third transform
      this->mTransform3D.init(spSetup3D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem3D = this->mTransform3D.requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORMS, mem3D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM3D, mem3D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::initCommunicator(SharedResolution spRes)
   {
      // initialise the communicator
      this->mCommunicator.init(spRes->spFwdSetup(Dimensions::Transform::TRA1D), spRes->spBwdSetup(Dimensions::Transform::TRA1D), spRes->spFwdSetup(Dimensions::Transform::TRA2D), spRes->spBwdSetup(Dimensions::Transform::TRA2D), spRes->spFwdSetup(Dimensions::Transform::TRA3D), spRes->spBwdSetup(Dimensions::Transform::TRA3D));
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::requiredOptions(std::set<NonDimensional::Id>& list) const
   {
      // Get list from 2 other dimensions
      Transform2DCoordinator<T1D,T2D, TCommunicator>::requiredOptions(list);

      this->mTransform3D.requiredOptions(list, Dimensions::Transform::TRA3D);
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::setOptions(const std::map<NonDimensional::Id, MHDFloat>& options)
   {
      // Get list from 2 other dimensions
      Transform2DCoordinator<T1D,T2D, TCommunicator>::setOptions(options);

      this->mTransform3D.setOptions(options, Dimensions::Transform::TRA3D);
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> std::vector<Array> Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::mesh()
   {
      std::vector<Array>   mesh;

      // Add grid of first dimension
      mesh.push_back(this->transform1D().meshGrid());

      // Add grid of second dimension
      mesh.push_back(this->transform2D().meshGrid());

      // Add grid of third dimension
      mesh.push_back(this->transform3D().meshGrid());

      return mesh;
   }

}

#endif // TRANSFORM3DCOORDINATOR_HPP

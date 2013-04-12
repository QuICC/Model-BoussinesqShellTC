/** \file Transform2DCoordinator.hpp
 *  \brief Implementation of a transform coordinator in 2D
 */

#ifndef TRANSFORM2DCOORDINATOR_HPP
#define TRANSFORM2DCOORDINATOR_HPP

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
#include "TransformCoordinators/Transform1DCoordinator.hpp"
#include "Variables/VariableRequirement.hpp"

namespace GeoMHDiSCC {

   /**
    * 
    * @brief Implementation of a transform coordinator in 2D
    *
    * The transform coordinator overlooks the whole transform process. Providing storage, communication
    * and the actual transforms
    */ 
   template <typename T1D, typename T2D, typename TCommunicator> class Transform2DCoordinator: public Transform1DCoordinator<T1D, TCommunicator>
   {
      public:
         /// Typedef for the first dimension transform
         typedef T1D Transform1DType;

         /// Typedef for the second dimension transform
         typedef T2D Transform2DType;

         /// Typedef for the communicator type
         typedef TCommunicator CommunicatorType;

         /**
          * @brief Very basic constructor
          */
         Transform2DCoordinator();

         /**
          * @brief Destructor
          */
         virtual ~Transform2DCoordinator();

         /**
          * @brief Initialise transforms of the 2D coordinator
          *
          * @param spRes   Resolution information object
          * @param varInfo Variables information
          */
         void initTransforms(SharedResolution spRes, const VariableRequirement& varInfo);

         /**
          * @brief Initialise the data communicator
          *
          * @param spRes   Resolution information object
          */
         void initCommunicator(SharedResolution spRes);

         /**
          * @brief Get the transform for the second dimension
          */
         T2D&   transform2D();

         /**
          * @brief Get grid array(s) of the mesh
          */
         std::vector<Array>  mesh();
         
      protected:

      private:
         /**
          * @brief The second dimension transform
          */
         T2D  mTransform2D;

         /**
          * @brief Initialise the transforms
          *
          * @param spSetup2D Shared setup object for 2D transform
          */
         void initTransform(typename T2D::SharedSetupType spSetup2D);
   };

   template <typename T1D, typename T2D, typename TCommunicator> inline T2D&  Transform2DCoordinator<T1D, T2D, TCommunicator>::transform2D()
   {
      return this->mTransform2D;
   }

   template <typename T1D, typename T2D, typename TCommunicator> Transform2DCoordinator<T1D, T2D, TCommunicator>::Transform2DCoordinator()
      : Transform1DCoordinator<T1D, TCommunicator>()
   {
   }

   template <typename T1D, typename T2D, typename TCommunicator> Transform2DCoordinator<T1D, T2D, TCommunicator>::~Transform2DCoordinator()
   {
   }

   template <typename T1D, typename T2D, typename TCommunicator> void Transform2DCoordinator<T1D, T2D, TCommunicator>::initTransforms(SharedResolution spRes, const VariableRequirement& varInfo)
   {
      // initialise the other dimension
      Transform1DCoordinator<T1D, TCommunicator>::initTransforms(spRes, varInfo);

      // Initialise the transforms
      this->initTransform(spRes->spTransformSetup(Dimensions::Transform::TRA2D));
   }

   template <typename T1D, typename T2D, typename TCommunicator> void Transform2DCoordinator<T1D, T2D, TCommunicator>::initTransform(typename T2D::SharedSetupType spSetup2D)
   {
      // initialise the second
      this->mTransform2D.init(spSetup2D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem2D = this->mTransform2D.requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORMS, mem2D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM2D, mem2D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename T1D, typename T2D, typename TCommunicator> void Transform2DCoordinator<T1D, T2D, TCommunicator>::initCommunicator(SharedResolution spRes)
   {
      // initialise the communicator
      this->mCommunicator.init(spRes->spFwdSetup(Dimensions::Transform::TRA1D), spRes->spBwdSetup(Dimensions::Transform::TRA1D), spRes->spFwdSetup(Dimensions::Transform::TRA2D), spRes->spBwdSetup(Dimensions::Transform::TRA2D));
   }

   template <typename T1D, typename T2D, typename TCommunicator> std::vector<Array> Transform2DCoordinator<T1D, T2D, TCommunicator>::mesh()
   {
      std::vector<Array>   mesh;

      // Add grid of first dimension
      mesh.push_back(this->transform1D().meshGrid());

      // Add grid of second dimension
      mesh.push_back(this->transform2D().meshGrid());

      return mesh;
   }

}

#endif // TRANSFORM2DCOORDINATOR_HPP
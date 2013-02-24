/** \file Transform3DCoordinator.hpp
 *  \brief Implementation of a transform coordinator in 3D
 */

#ifndef TRANSFORM3DCOORDINATOR_HPP
#define TRANSFORM3DCOORDINATOR_HPP

// Debug includes
//
#include "Debug/PrepMacros/StorageProfilerMacro.h"

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Transforms/Transform2DCoordinator.hpp"

namespace GeoMHDiSCC {

   /**
    * 
    * @brief Implementation of a transform coordinator in 3D
    *
    * The transform coordinator overlooks the whole transform process. Providing storage, communication
    * and the actual transforms
    */ 
   template <typename T1D, typename T2D, typename T3D, typename TCommunicator = CommunicatorMacro3D> class Transform3DCoordinator: public Transform2DCoordinator<T1D, T2D, TCommunicator>
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
          * @param setup1D Setup object for 1D transform
          * @param setup2D Setup object for 2D transform
          * @param setup3D Setup object for 3D transform
          * @param varInfo Variables information
          */
         void initTransforms(const typename T1D::SetupType& setup1D, const typename T2D::SetupType& setup2D, const typename T3D::SetupType& setup3D, const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Initialise the data communicator
          *
          * @param spRes   Resolution information object
          */
         void initCommunicator(SharedResolution spRes);

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
          * @param setup3D Setup object for 3D transform
          */
         void initTransform(const typename T3D::SetupType& setup3D);
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

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::initTransforms(const typename T1D::SetupType& setup1D, const typename T2D::SetupType& setup2D, const typename T3D::SetupType& setup3D, const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {
      // initialise 2 other dimensions
      Transform2DCoordinator<T1D,T2D, TCommunicator>::initTransforms(setup1D, setup2D, varInfo);

      // Initialise the transforms
      this->initTransform(setup3D);
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::initTransform(const typename T3D::SetupType& setup3D)
   {
      // initialise the third transform
      this->mTransform3D.init(setup3D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem3D = this->mTransform3D.requiredStorage();
         StorageProfilerMacro_update(StorageProfiler::TRANSFORMS, mem3D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            DetailedStorageProfilerMacro_update(StorageProfiler::TRANSFORM3D, mem3D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename T1D, typename T2D, typename T3D, typename TCommunicator> void Transform3DCoordinator<T1D, T2D, T3D, TCommunicator>::initCommunicator(SharedResolution spRes)
   {
      // initialise the communicator
      this->mCommunicator.init(*spRes->spFwdSetup(0), *spRes->spBwdSetup(0), *spRes->spFwdSetup(1), *spRes->spBwdSetup(1), *spRes->spFwdSetup(2), *spRes->spBwdSetup(2));
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

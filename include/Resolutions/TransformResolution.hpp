/** \file TransformResolution.hpp
 *  \brief Definition of the required resolution information for a transform
 */

#ifndef TRANSFORMRESOLUTION_HPP
#define TRANSFORMRESOLUTION_HPP

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
#include "ScalarFields/ScalarFieldSetup.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Definition of the required resolution information for a transform
    */
   class TransformResolution
   {
      public:
         /**
          * @brief Empty constructor
          *
          * @param fwd Set of indexes for 1D (forward direction)
          * @param bwd Set of indexes for 1D (backward direction)
          * @param idx2D Set of indexes for 2D
          * @param idx3D Set of indexes for 3D
          */
         TransformResolution(const std::vector<ArrayI>& fwd, const std::vector<ArrayI>& bwd, const std::vector<ArrayI>& idx2D, const ArrayI& idx3D);

         /**
          * @brief Empty Destructor
          */
         virtual ~TransformResolution();

         /**
          * @brief Get first dimension (forward direction)
          *
          * @param j Index of the second dimension
          * @param k Index of the third dimension
          */
         int dimFwd(const int j = 0, const int k = 0) const;

         /**
          * @brief Get first dimension (backward direction)
          *
          * @param j Index of the second dimension
          * @param k Index of the third dimension
          */
         int dimBwd(const int j = 0, const int k = 0) const;

         /**
          * @brief Get second dimension
          *
          * @param k Index of the third dimension 
          */
         int dim2D(const int k = 0) const;

         /**
          * @brief Get third dimension
          */
         int dim3D() const;

         /**
          * @brief Get index mapping to full resolution for 1D (forward direction)
          *
          * @param i Index of the first dimension
          * @param j Index of the second dimension
          * @param k Index of the third dimension
          */
         int idxFwd(const int i, const int j = 0, const int k = 0) const;

         /**
          * @brief Get index mapping to full resolution for 1D (backward direction)
          *
          * @param i Index of the first dimension
          * @param j Index of the second dimension
          * @param k Index of the third dimension
          */
         int idxBwd(const int i, const int j = 0, const int k = 0) const;

         /**
          * @brief Get index mapping to full resolution for 2D
          *
          * @param j Index of the second dimension
          * @param k Index of the third dimension 
          */
         int idx2D(const int j, const int k = 0) const;

         /**
          * @brief Get full set of index mapping to full resolution for 3D
          */
         const ArrayI& idx3D() const;

         /**
          * @brief Get index mapping to full resolution for 3D
          *
          * @param k Index of the third dimension 
          */
         int idx3D(const int k) const;

         /**
          * @brief Get the forward scalar field setup
          */
         Datatypes::SharedScalarFieldSetup spFwdSetup() const;

         /**
          * @brief Get the backward scalar field setup
          */
         Datatypes::SharedScalarFieldSetup spBwdSetup() const;

      protected:

      private:
         /**
          * @brief Set of indexes describing the first dimensions (forward direction)
          */
         std::vector<ArrayI>   mFwd;

         /**
          * @brief Set of indexes describing the first dimensions (backward direction)
          */
         std::vector<ArrayI>   mBwd;

         /**
          * @brief Set of indexes describing the second dimensions
          */
         std::vector<ArrayI>   mIdx2D;

         /**
          * @brief Set of indexes describing the third dimensions
          */
         ArrayI   mIdx3D;
   };

   /// Typedef for a shared pointer to a TransformResolution object
   typedef SharedPtrMacro<TransformResolution>   SharedTransformResolution;

   /// Typedef for a shared pointer to a const TransformResolution object
   typedef SharedPtrMacro<const TransformResolution>   SharedCTransformResolution;

}

#endif // TRANSFORMRESOLUTION_HPP

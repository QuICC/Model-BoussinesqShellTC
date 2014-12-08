/** 
 * @file CylindricalTransformSteps.cpp
 * @brief Source of the implementation of the cylindrical transform steps
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/CylindricalTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<TransformBranch>  backwardScalar()
   {
      std::vector<TransformBranch> transform;

      transform.push_back(FieldComponents::Spectral::SCALAR, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR);

      return transform;
   }

   std::vector<TransformBranch>  backwardGradient()
   {
      std::vector<TransformBranch> transform;

      transform.push_back(FieldComponents::Spectral::SCALAR, TransformBranch::Proj1DType::DIFF, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::SCALAR, TransformBranch::Proj1DType::DIVR, TransformBranch::Proj2DType::DIFF, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::SCALAR, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT);

      return transform;
   }

   std::vector<TransformBranch>  backwardVector()
   {
      std::vector<TransformBranch> transform;

      transform.push_back(FieldComponents::Spectral::ONE, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::VECTOR);
      transform.push_back(FieldComponents::Spectral::TWO, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::VECTOR);
      transform.push_back(FieldComponents::Spectral::THREE, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::VECTOR);

      return transform;
   }

   std::vector<TransformBranch>  backwardVGradient()
   {
      std::vector<TransformBranch> transform;

      transform.push_back(FieldComponents::Spectral::ONE, TransformBranch::Proj1DType::DIFF, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::ONE, TransformBranch::Proj1DType::DIVR, TransformBranch::Proj2DType::DIFF, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::ONE, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT);

      transform.push_back(FieldComponents::Spectral::TWO, TransformBranch::Proj1DType::DIFF, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::TWO, TransformBranch::Proj1DType::DIVR, TransformBranch::Proj2DType::DIFF, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::TWO, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT);

      transform.push_back(FieldComponents::Spectral::THREE, TransformBranch::Proj1DType::DIFF, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::THREE, TransformBranch::Proj1DType::DIVR, TransformBranch::Proj2DType::DIFF, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT);
      transform.push_back(FieldComponents::Spectral::THREE, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT);

      return transform;
   }

   std::vector<TransformBranch>  backwardCurl()
   {
      std::vector<TransformBranch> transform;

      transform.push_back(FieldComponents::Spectral::ONE, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::DIFF, FieldComponents::Physical::TWO, FieldType::CURL);
      transform.push_back(FieldComponents::Spectral::ONE, TransformBranch::Proj1DType::DIVR, TransformBranch::Proj2DType::DIFF, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL);

      transform.push_back(FieldComponents::Spectral::TWO, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::DIFF, FieldComponents::Physical::ONE, FieldType::CURL);
      transform.push_back(FieldComponents::Spectral::TWO, TransformBranch::Proj1DType::DIVRDIFFR, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL);

      transform.push_back(FieldComponents::Spectral::THREE, TransformBranch::Proj1DType::DIVR, TransformBranch::Proj2DType::DIFF, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::CURL);
      transform.push_back(FieldComponents::Spectral::THREE, TransformBranch::Proj1DType::DIFF, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::CURL);

      return transform;
   }

   std::vector<TransformBranch>  backwardDivergence()
   {
      std::vector<TransformBranch> transform;

      transform.push_back(FieldComponents::Spectral::ONE, TransformBranch::Proj1DType::DIVRDIFFR, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE);

      transform.push_back(FieldComponents::Spectral::TWO, TransformBranch::Proj1DType::DIVR, TransformBranch::Proj2DType::DIFF, TransformBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE);

      transform.push_back(FieldComponents::Spectral::THREE, TransformBranch::Proj1DType::PROJ, TransformBranch::Proj2DType::PROJ, TransformBranch::Proj3DType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE);

      return transform;
   }

}
}
}

/** 
 * @file Cartesian2DTransformSteps.cpp
 * @brief Source of the implementation of the 2D cartesian transform steps
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
#include "TransformConfigurators/Cartesian2DTransformSteps.hpp"

// Project includes
//

namespace QuICC {

namespace Transform {

namespace TransformSteps {

   std::vector<TransformPath>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<TransformPath> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      if(isNL)
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::ADD);
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 2);
      std::vector<TransformPath> transform;

      if(isNL)
      {
         transform.push_back(TransformPath(FieldComponents::Physical::ONE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::ONE, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::TWO, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::TWO, Arithmetics::ADD);
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::ONE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::ONE, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::TWO, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::TWO, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      return transform;
   }

   std::vector<TransformPath>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardGradient2(const std::map<std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>,bool>& req)
   {
      std::vector<TransformPath> transform;
      std::pair<FieldComponents::Physical::Id,FieldComponents::Physical::Id>  pairId;

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::ONE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::DIFF2);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(ProjectorNDType::DIFF, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF2, pairId, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      throw Exception("Curl cannot be computed in 2D!");

      std::vector<TransformPath> transform;

      return transform;
   }

   std::vector<TransformPath>  backwardDivergence()
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIFF);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      return transform;
   }

}
}
}

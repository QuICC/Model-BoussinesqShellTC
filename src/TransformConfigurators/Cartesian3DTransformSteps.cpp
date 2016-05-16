/** 
 * @file Cartesian3DTransformSteps.cpp
 * @brief Source of the implementation of the 3D cartesian transform steps
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
#include "TransformConfigurators/Cartesian3DTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

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
         transform.back().addEdge(IntgeratorNDType::INTG, Integrator2DType::INTG);
         transform.back().addEdge(Integratro2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::ADD);
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::SCALAR, FieldType::SCALAR));
         transform.back().addEdge(IntgeratorNDType::INTG, Integrator2DType::INTG);
         transform.back().addEdge(Integratro2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, scalId, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 3);
      std::vector<TransformPath> transform;

      if(isNL)
      {
         transform.push_back(TransformPath(FieldComponents::Physical::ONE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::ONE, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::TWO, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::TWO, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::THREE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::THREE, Arithmetics::ADD);
      } else
      {
         transform.push_back(TransformPath(FieldComponents::Physical::ONE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::ONE, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::TWO, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::TWO, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Physical::THREE, FieldType::VECTOR));
         transform.back().addEdge(IntegratorNDType::INTG);
         transform.back().addEdge(Integrator2DType::INTG);
         transform.back().addEdge(Integrator1DType::INTG, FieldComponents::Spectral::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(Projector2DType::PROJ);
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
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THREE, Arithmetics::ADD);
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
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::ONE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::TWO);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF2);
         transform.back().addEdge(ProjectorNDType::PROJ, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::TWO,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, pairId, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::DIFF, pairId, Arithmetics::ADD);
      }

      pairId = std::make_pair(FieldComponents::Physical::THREE,FieldComponents::Physical::THREE);
      if(req.find(pairId)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::SCALAR, FieldType::GRADIENT2));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
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
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::VECTOR));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THREE, Arithmetics::ADD);
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
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(id, FieldType::GRADIENT));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<TransformPath> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::ONE, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::ONE, Arithmetics::ADD);
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::TWO, Arithmetics::ADD);

         transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::TWO, Arithmetics::SUB);
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::CURL));
         transform.back().addEdge(Projector1DType::PROJ);
         transform.back().addEdge(Projector2DType::DIFF);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THREE, Arithmetics::SUB);

         transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::CURL));
         transform.back().addEdge(Projector1DType::DIFF);
         transform.back().addEdge(Projector2DType::PROJ);
         transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::THREE, Arithmetics::ADD);
      }

      return transform;
   }

   std::vector<TransformPath>  backwardDivergence()
   {
      std::vector<TransformPath> transform;

      transform.push_back(TransformPath(FieldComponents::Spectral::ONE, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::DIFF);
      transform.back().addEdge(Projector2DType::PROJ);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::TWO, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(Projector2DType::DIFF);
      transform.back().addEdge(ProjectorNDType::PROJ, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      transform.push_back(TransformPath(FieldComponents::Spectral::THREE, FieldType::DIVERGENCE));
      transform.back().addEdge(Projector1DType::PROJ);
      transform.back().addEdge(Projector2DType::PROJ);
      transform.back().addEdge(ProjectorNDType::DIFF, FieldComponents::Physical::SCALAR, Arithmetics::ADD);

      return transform;
   }

}
}
}

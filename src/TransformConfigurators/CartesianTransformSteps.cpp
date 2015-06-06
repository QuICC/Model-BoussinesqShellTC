/** 
 * @file CartesianTransformSteps.cpp
 * @brief Source of the implementation of the cartesian transform steps
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
#include "TransformConfigurators/CartesianTransformSteps.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<IntegratorBranch>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<IntegratorBranch> transform;

      FieldComponents::Spectral::Id scalId = components.at(0).first;

      if(isNL)
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, scalId, FieldType::SCALAR));
      } else
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, scalId, FieldType::SCALAR));
      }

      return transform;
   }

   std::vector<IntegratorBranch>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 3);
      std::vector<IntegratorBranch> transform;

      if(isNL)
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::ONE, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::ONE, FieldType::VECTOR));

         transform.push_back(IntegratorBranch(FieldComponents::Physical::TWO, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::TWO, FieldType::VECTOR));

         transform.push_back(IntegratorBranch(FieldComponents::Physical::THREE, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::THREE, FieldType::VECTOR));
      } else
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::ONE, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::ONE, FieldType::VECTOR));

         transform.push_back(IntegratorBranch(FieldComponents::Physical::TWO, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::TWO, FieldType::VECTOR));

         transform.push_back(IntegratorBranch(FieldComponents::Physical::THREE, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::THREE, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR));

      return transform;
   }

   std::vector<ProjectorBranch>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THREE, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::ONE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::ONE, FieldType::CURL));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::ONE, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::TWO)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::TWO, FieldType::CURL));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::TWO, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::THREE)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL));

         transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THREE, FieldType::CURL));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardDivergence()
   {
      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::ONE, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::TWO, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::THREE, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

}
}
}
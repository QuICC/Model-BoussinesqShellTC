/** 
 * @file CylinderTransformSteps.cpp
 * @brief Source of the implementation of the phyiscal <-> spectral transform steps in a whole cylinder
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
#include "TransformConfigurators/CylinderTransformSteps.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Transform {

namespace TransformSteps {

   std::vector<IntegratorBranch>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<IntegratorBranch> transform;

      if(isNL)
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      } else
      {
         transform.push_back(IntegratorBranch(FieldComponents::Physical::SCALAR, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      }

      return transform;
   }

   std::vector<IntegratorBranch>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 2);
      std::vector<IntegratorBranch> transform;
      FieldComponents::Spectral::Id radialId = components.at(0).first;
      int radialFlag = components.at(0).second;
      FieldComponents::Spectral::Id thetaId = components.at(1).first;
      int thetaFlag = components.at(1).second;
      FieldComponents::Spectral::Id zId = components.at(1).first;
      int zFlag = components.at(1).second;

      if(isNL)
      {
         if(radialFlag == 0)
         {
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, radialId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(thetaFlag == 0)
         {
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, thetaId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(zFlag == 0)
         {
            transform.push_back(IntegratorBranch(FieldComponents::Physical::Z, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, zId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }
      } else
      {
         if(radialFlag == 0)
         {
            transform.push_back(IntegratorBranch(FieldComponents::Physical::R, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, radialId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(thetaFlag == 0)
         {
            transform.push_back(IntegratorBranch(FieldComponents::Physical::THETA, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, thetaId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(zFlag == 0)
         {
            transform.push_back(IntegratorBranch(FieldComponents::Physical::Z, IntegratorBranch::Intg3DType::INTG, IntegratorBranch::Intg2DType::INTG, IntegratorBranch::Intg1DType::INTG, zId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }
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

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::SCALAR, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::Z, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::Z, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::Z, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch> transform;

      if(id == FieldComponents::Spectral::R || id == FieldComponents::Spectral::THETA)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIFFDIVR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIVR2, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
         }
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch(id, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::Z, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      throw Exception("Implementation of the curl is not correct!");

      std::vector<ProjectorBranch> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::R, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::Z, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::R, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::Z, ProjectorBranch::Proj1DType::DIFF, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::THETA, FieldType::CURL));
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::Z, FieldType::CURL));
         transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::Z, FieldType::CURL));
      }

      return transform;
   }

   std::vector<ProjectorBranch>  backwardDivergence()
   {
      throw Exception("Implementation of the divergence is not correct!");

      std::vector<ProjectorBranch> transform;

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::R, ProjectorBranch::Proj1DType::DIVRDIFFR, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::THETA, ProjectorBranch::Proj1DType::DIVR, ProjectorBranch::Proj2DType::DIFF, ProjectorBranch::Proj3DType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch(FieldComponents::Spectral::Z, ProjectorBranch::Proj1DType::PROJ, ProjectorBranch::Proj2DType::PROJ, ProjectorBranch::Proj3DType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

}
}
}

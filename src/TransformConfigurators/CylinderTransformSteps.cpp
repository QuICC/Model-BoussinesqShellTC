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

   std::vector<IntegratorBranch3D>  forwardScalar(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 1);
      std::vector<IntegratorBranch3D> transform;

      if(isNL)
      {
         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      } else
      {
         transform.push_back(IntegratorBranch3D(FieldComponents::Physical::SCALAR, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, FieldComponents::Spectral::SCALAR, FieldType::SCALAR));
      }

      return transform;
   }

   std::vector<IntegratorBranch3D>  forwardVector(const std::vector<std::pair<FieldComponents::Spectral::Id,int> >& components, const bool isNL)
   {
      assert(components.size() == 3);
      std::vector<IntegratorBranch3D> transform;
      FieldComponents::Spectral::Id radialId = components.at(0).first;
      int radialFlag = components.at(0).second;
      FieldComponents::Spectral::Id thetaId = components.at(1).first;
      int thetaFlag = components.at(1).second;
      FieldComponents::Spectral::Id zId = components.at(2).first;
      int zFlag = components.at(2).second;

      if(isNL)
      {
         if(radialFlag == 0)
         {
            transform.push_back(IntegratorBranch3D(FieldComponents::Physical::R, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, radialId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(thetaFlag == 0)
         {
            transform.push_back(IntegratorBranch3D(FieldComponents::Physical::THETA, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, thetaId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(zFlag == 0)
         {
            transform.push_back(IntegratorBranch3D(FieldComponents::Physical::Z, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, zId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }
      } else
      {
         if(radialFlag == 0)
         {
            transform.push_back(IntegratorBranch3D(FieldComponents::Physical::R, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, radialId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(thetaFlag == 0)
         {
            transform.push_back(IntegratorBranch3D(FieldComponents::Physical::THETA, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, thetaId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }

         if(zFlag == 0)
         {
            transform.push_back(IntegratorBranch3D(FieldComponents::Physical::Z, IntgPhysType::INTG, IntgPartType::INTG, IntgSpecType::INTG, zId, FieldType::VECTOR));
         } else
         {
            throw Exception("Requested an unknown vector forward transform");
         }
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardScalar(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::SCALAR));

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardGradient(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::DIVR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::SCALAR, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::Z, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardVector(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::R, ProjSpecType::DIVR, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THETA, ProjSpecType::DIVR, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::VECTOR));
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::Z, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::Z, FieldType::VECTOR));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardVGradient(FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(id == FieldComponents::Spectral::R || id == FieldComponents::Spectral::THETA)
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(ProjectorBranch3D(id, ProjSpecType::DIFFDIVR, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(ProjectorBranch3D(id, ProjSpecType::DIVR2, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
         }
      } else
      {
         if(req.find(FieldComponents::Physical::R)->second)
         {
            transform.push_back(ProjectorBranch3D(id, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::GRADIENT));
         }

         if(req.find(FieldComponents::Physical::THETA)->second)
         {
            transform.push_back(ProjectorBranch3D(id, ProjSpecType::DIVR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::GRADIENT));
         }
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch3D(id, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::Z, FieldType::GRADIENT));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardCurl(const std::map<FieldComponents::Physical::Id,bool>& req)
   {
      std::vector<ProjectorBranch3D> transform;

      if(req.find(FieldComponents::Physical::R)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THETA, ProjSpecType::DIVR, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::R, FieldType::CURL, Arithmetics::SETNEG));

         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::Z, ProjSpecType::DIVR, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::R, FieldType::CURL, Arithmetics::ADD));
      }

      if(req.find(FieldComponents::Physical::THETA)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::R, ProjSpecType::DIVR, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::SET));

         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::Z, ProjSpecType::DIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::THETA, FieldType::CURL, Arithmetics::SUB));
      }

      if(req.find(FieldComponents::Physical::Z)->second)
      {
         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::R, ProjSpecType::DIVR2, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::Z, FieldType::CURL, Arithmetics::SETNEG));

         transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THETA, ProjSpecType::DIVRDIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::Z, FieldType::CURL, Arithmetics::ADD));
      }

      return transform;
   }

   std::vector<ProjectorBranch3D>  backwardDivergence()
   {
      std::vector<ProjectorBranch3D> transform;

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::R, ProjSpecType::DIVRDIFF, ProjPartType::PROJ, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::THETA, ProjSpecType::DIVR2, ProjPartType::DIFF, ProjPhysType::PROJ, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      transform.push_back(ProjectorBranch3D(FieldComponents::Spectral::Z, ProjSpecType::PROJ, ProjPartType::PROJ, ProjPhysType::DIFF, FieldComponents::Physical::SCALAR, FieldType::DIVERGENCE));

      return transform;
   }

}
}
}
